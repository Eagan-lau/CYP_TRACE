import unittest
from build_multiaxis_splits import chemical_components, fragment_keys, make_block, fold_map


class MultiAxisTest(unittest.TestCase):
    def test_both_sides_and_coproduct_bridges_preserved(self):
        rows={'a':{'substrates':['CC'],'products':['CCO','C=O']},
              'b':{'substrates':['CCC'],'products':['CCCO','C=O']},
              'c':{'substrates':['c1ccccc1'],'products':['Oc1ccccc1']}}
        keys,groups,profile=chemical_components(rows)
        self.assertEqual(sorted(sorted(v) for v in groups.values()),[['a','b'],['c']])
        self.assertIn('acyclic_connectivity:C=O',keys['a'])
        self.assertEqual(profile['acyclic_connectivity:C=O'],{'a','b'})

    def test_stereo_collapses_only_in_scaffold_key_and_acyclic_not_empty(self):
        self.assertEqual(fragment_keys('C[C@H](O)CC'),fragment_keys('C[C@@H](O)CC'))
        self.assertNotEqual(fragment_keys('CC'),fragment_keys('CCC'))
        self.assertNotEqual(fragment_keys('CC(=O)O'),fragment_keys('CC(=O)[O-]'))

    def test_double_cold_and_publication_scope(self):
        rows=[{'sequence_sha256':'p0','reaction_key':'c0'},
              {'sequence_sha256':'p0','reaction_key':'c1'},
              {'sequence_sha256':'p1','reaction_key':'c1'},
              {'sequence_sha256':'p2','reaction_key':'c1'}]
        units=[{'A'},{'B'},{'B'},{'C'}]
        block=make_block('double_cold',list(range(4)),rows,units,
            {'p0':0,'p1':1,'p2':1},{'c0':0,'c1':1},0,0,
            {'c0':['key0'],'c1':['key1']},{'p0':'g0','p1':'g1','p2':'g2'},
            {'c0':'h0','c1':'h1'},2)
        self.assertEqual(block['test_edge_indices'],[0])
        self.assertEqual(block['publication_scope_edge_indices'],[0,1])
        self.assertEqual(block['train_edge_indices'],[3])
        self.assertEqual(block['publication_purged_train_edges'],1)

    def test_empty_training_not_replaced(self):
        rows=[{'sequence_sha256':'p0','reaction_key':'c0'},{'sequence_sha256':'p1','reaction_key':'c1'}]
        block=make_block('double_cold',[0,1],rows,[{'A'},{'A'}],{'p0':0,'p1':1},
            {'c0':0,'c1':1},0,0,{'c0':['k0'],'c1':['k1']},
            {'p0':'g0','p1':'g1'},{'c0':'h0','c1':'h1'},2)
        self.assertFalse(block['evaluable'])
        self.assertEqual(block['train_edge_indices'],[])

    def test_fold_map_never_splits_existing_group(self):
        mapping,k=fold_map({'g0':['a','b','c'],'g1':['d'],'g2':['e']},{'a','b','c','d','e'},3)
        self.assertEqual(k,3)
        self.assertEqual(len({mapping[s] for s in ['a','b','c']}),1)


if __name__=='__main__': unittest.main()
