"""Editable Word equations for exact, audited Markdown math tokens.

Changes presentation only; the authoritative Markdown expressions are unchanged.
Unrecognised code tokens continue to use the original code renderer.
"""
import re
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH


def run(text, roman=False):
    r = OxmlElement("m:r")
    if roman:
        pr = OxmlElement("m:rPr")
        # Literal normal text prevents LibreOffice's math parser treating words
        # such as "end" and the hyphens in an index as equation commands.
        for tag in ("lit", "nor"):
            flag=OxmlElement("m:"+tag)
            flag.set(qn("m:val"), "1")
            pr.append(flag)
        sty = OxmlElement("m:sty")
        sty.set(qn("m:val"), "p")
        pr.append(sty)
        r.append(pr)
    wp = OxmlElement("w:rPr")
    font = OxmlElement("w:rFonts")
    for attr in ("ascii", "hAnsi"):font.set(qn("w:" + attr), "Cambria Math")
    wp.append(font)
    sz = OxmlElement("w:sz")
    sz.set(qn("w:val"), "20")
    wp.append(sz)
    r.append(wp)
    t = OxmlElement("m:t")
    t.set(qn("xml:space"), "preserve")
    t.text = text
    r.append(t)
    return r


def field(tag, *nodes):
    e=OxmlElement("m:"+tag)
    e.extend(nodes)
    return e


def sub(symbol, index, roman=False):
    return field("sSub", field("e", run(symbol, roman)), field("sub", run(index, True)))


def power(node, exp):
    return field("sSup", field("e", node), field("sup", run(exp)))


def fraction(num, den):
    return field("f", field("num", *num), field("den", *den))


def weighted_sum(symbols):
    return [sub("∑", "i", True), *[sub(s, "i") for s in symbols]]


def identity(end="end-to-end"):
    return [sub("MRR", end, True), run(" = "), sub("C", "w"), run(" × "), sub("MRR", "covered,w", True)]


def expression(token):
    if token=="MRR_end-to-end = C_w × MRR_covered,w":return identity()
    if token=="score(p,c) = log P_chem(c) + lambda × Delta(p,c)":
        return [run("score", True),run("(p,c) = "),run("log ",True),sub("P","chem"),run("(c) + λ × Δ(p,c)")]
    if token=="Delta_int(p,r) = X_p W Z_r^T":
        return [sub("Δ","int"),run("(p,r) = "),sub("X","p"),run(" W "),power(sub("Z","r"),"T")]
    if token in ("X","Z"):return [run(token)]
    if token in ("s_i","t_j"):return [sub(*token.split("_"))]
    if token=="alpha":return [run("α")]
    if token=="s_i^2 t_j^2 / (s_i^2 t_j^2 + alpha)":
        pair=lambda:[power(sub("s","i"),"2"),power(sub("t","j"),"2")]
        return [fraction(pair(),pair()+[run(" + α")])]
    return None


DECOMPOSITION="C_w = sum_i(w_i a_i) / sum_i(w_i);  MRR_covered,w = sum_i(w_i a_i r_i) / sum_i(w_i a_i);  MRR_end-to-end,w = C_w × MRR_covered,w."


def append_equation(paragraph, nodes):
    paragraph._p.append(field("oMath",*nodes))


def install(base):
    original_inline=base.add_inline
    original_body=base.add_body_paragraph

    def add_inline(paragraph,text,base_size=10.0):
        position=0
        for match in re.finditer(r"`([^`]+)`",text):
            nodes=expression(match.group(1))
            if nodes is None:continue
            original_inline(paragraph,text[position:match.start()],base_size)
            append_equation(paragraph,nodes)
            position=match.end()
        original_inline(paragraph,text[position:],base_size)

    def add_body(document,text):
        if text=="`"+DECOMPOSITION+"`":
            p=document.add_paragraph(style="Normal")
            p.alignment=WD_ALIGN_PARAGRAPH.CENTER
            equations=[
                [sub("C","w"),run(" = "),fraction(weighted_sum("wa"),weighted_sum("w"))],
                [sub("MRR","covered,w",True),run(" = "),fraction(weighted_sum("war"),weighted_sum("wa"))],
                identity("end-to-end,w")+[run(".")],
            ]
            for i,nodes in enumerate(equations):
                if i:p.add_run().add_break()
                append_equation(p,nodes)
            p.paragraph_format.keep_together=True
            return p
        return original_body(document,text)

    base.add_inline=add_inline
    base.add_body_paragraph=add_body
