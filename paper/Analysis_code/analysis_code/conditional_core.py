"""Train-only linear kernel residual plus convex, group-weighted OOF selection."""
from collections import Counter
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar, brentq
from scipy.special import logsumexp


def normalized_features(features):
    x=np.asarray(features,dtype=np.float64)
    if x.ndim!=2 or not np.all(np.isfinite(x)): raise ValueError('Invalid protein features')
    return x/np.maximum(np.linalg.norm(x,axis=1,keepdims=True),1e-12)


def fit_residual(P,Y,K):
    training=np.flatnonzero(Y.sum(1))
    if not len(training): raise ValueError('Empty training pool')
    balanced=Y[training].astype(np.float64)/Y[training].sum(1,keepdims=True)
    target=balanced@np.asarray(K,dtype=np.float64); prior=target.mean(0)
    probability=np.maximum(prior/prior.sum(),1e-12); probability/=probability.sum()
    residual=target-prior; rms=float(np.sqrt(np.mean(residual**2)))
    x=P[training]; center=x.mean(0); xc=x-center
    pk=xc@xc.T; kscale=float(np.trace(pk)/len(training))
    result={'training':training,'center':center,'training_features_centered':xc,'kernel_scale':kscale,
        'prior_log':np.log(probability),'residual_rms':rms,'alpha':None,'alpha_bound_adjacent':False,
        'zero_residual':True,'zero_reason':None,'coefficients':np.zeros((len(training),K.shape[1]))}
    if len(training)<3 or rms<1e-12 or kscale<1e-12:
        result['zero_reason']='insufficient_training_variation'; return result
    pk/=kscale
    eigen,u=eigh(pk,check_finite=False); eigen=np.maximum(eigen,0)
    z=u.T@residual; energy=np.sum(z*z,axis=1); n=len(training)
    def gcv(logalpha):
        alpha=np.exp(logalpha); shrink=alpha/(eigen+alpha)
        denom=max(1-(1+np.sum(1-shrink))/n,1e-12)
        return float(np.sum(energy*shrink**2)/n/residual.shape[1]/denom**2)
    zero=float(np.mean(residual**2)/(1-1/n)**2)
    opt=minimize_scalar(gcv,bounds=(-14.,14.),method='bounded',options={'xatol':1e-6})
    result.update({'gcv_zero':zero,'gcv_selected':float(opt.fun),'alpha_bound_adjacent':abs(float(opt.x))>13.99})
    if not opt.success or not np.isfinite(opt.fun) or opt.fun>=zero*(1-1e-10):
        result['zero_reason']='gcv_favors_zero_residual'; return result
    alpha=float(np.exp(opt.x)); result['alpha']=alpha
    result['coefficients']=u@(z/(eigen[:,None]+alpha))
    result['zero_residual']=False
    return result


def predict_delta(model,Pquery):
    if model['zero_residual']: return np.zeros((len(Pquery),len(model['prior_log'])))
    cross=(Pquery-model['center'])@model['training_features_centered'].T/model['kernel_scale']
    return (cross@model['coefficients'])/model['residual_rms']


def group_weights(records):
    panels=Counter((r['group'],r['query']) for r in records)
    queries=Counter(g for g,q in panels)
    return np.array([1/len(queries)/queries[r['group']]/panels[(r['group'],r['query'])] for r in records])


def select_lambda(logprior,delta,targets,records):
    if len(records)!=len(logprior) or len(records)!=len(delta): raise ValueError('Calibration rows mismatch')
    groups={r['group'] for r in records}
    if len(groups)<2: return {'lambda':0.,'reason':'fewer_than_two_inner_groups','inner_groups':len(groups)}
    weights=group_weights(records)
    true_delta=np.array([np.mean(delta[i,truth]) for i,truth in enumerate(targets)])
    true_prior=np.array([np.mean(logprior[i,truth]) for i,truth in enumerate(targets)])
    def loss(value):
        return float(weights@(logsumexp(logprior+value*delta,axis=1)-true_prior-value*true_delta))
    def derivative(value):
        score=logprior+value*delta
        probability=np.exp(score-logsumexp(score,axis=1,keepdims=True))
        return float(weights@(np.sum(probability*delta,axis=1)-true_delta))
    base=loss(0.); d0=derivative(0.)
    receipt={'inner_groups':len(groups),'inner_rows':len(records),'zero_loss':base,'derivative_at_zero':d0}
    if d0>=-1e-12: return {**receipt,'lambda':0.,'selected_loss':base,'reason':'zero_boundary_optimum'}
    upper=1.
    while upper<1024 and derivative(upper)<=0: upper*=2
    if derivative(upper)<=0:
        return {**receipt,'lambda':0.,'selected_loss':base,'reason':'no_finite_bracket_identified'}
    value=float(brentq(derivative,0,upper,xtol=1e-8))
    return {**receipt,'lambda':value,'selected_loss':loss(value),'reason':'convex_inner_optimum','bracket_upper':upper}


def model_receipt(model):
    return {key:model[key] for key in ['alpha','alpha_bound_adjacent','zero_residual','zero_reason','residual_rms','kernel_scale']}
