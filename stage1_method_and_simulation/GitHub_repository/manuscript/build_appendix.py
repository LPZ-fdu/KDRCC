"""Build the Chinese appendix with editable Office Math equations.

Requires Python 3, python-docx, lxml and Microsoft's MML2OMML.XSL.
Usage: python build_appendix.py --xsl PATH [--output PATH] [--template PATH]
"""
import argparse, copy, json, re
from pathlib import Path
from lxml import etree as E
from docx import Document
from docx.shared import Cm, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

MM='http://www.w3.org/1998/Math/MathML'
SYMBOLS=dict(zip('alpha beta gamma delta epsilon varepsilon eta theta kappa lambda mu pi rho sigma tau phi varphi chi psi omega Delta Phi Psi Omega Sigma'.split(),
                 'α β γ δ ε ε η θ κ λ μ π ρ σ τ φ φ χ ψ ω Δ Φ Ψ Ω Σ'.split()))
SYMBOLS.update({'sum':'∑','prod':'∏','int':'∫','partial':'∂','nabla':'∇','in':'∈','notin':'∉','ne':'≠','le':'≤','ge':'≥','to':'→','mapsto':'↦','otimes':'⊗','cap':'∩','varnothing':'∅','cdots':'⋯','times':'×','pm':'±','sim':'∼','subseteq':'⊆','setminus':'∖','cup':'∪','mid':'|','ldots':'…','dagger':'†','prime':'′','lceil':'⌈','rceil':'⌉'})
OPS={'exp','log','min','max','sin','cos','tanh','arg'}
SYMBOLS['ell']='ℓ'
SYMBOLS['xi']='ξ'
SYMBOLS['Pi']='Π'
SYMBOLS['nu']='ν'
SYMBOLS.update({'infty':'∞','Rightarrow':'⇒','equiv':'≡'})
def el(tag,text=None,**attrs):
    x=E.Element('{'+MM+'}'+tag,**attrs)
    if text is not None:x.text=text
    return x
class MathParser:
    def __init__(self,s):self.s=s;self.i=0
    def skip(self):
        while self.i<len(self.s) and self.s[self.i].isspace():self.i+=1
    def seq(self,end=None):
        a=el('mrow')
        while True:
            self.skip()
            if self.i>=len(self.s):
                if end:raise ValueError('Unclosed group: '+self.s)
                break
            if end and self.s[self.i]==end:self.i+=1;break
            if self.s[self.i]=='}':raise ValueError('Unexpected brace: '+self.s)
            base=self.atom();lo=hi=None
            while True:
                self.skip()
                if self.i>=len(self.s) or self.s[self.i] not in '_^':break
                c=self.s[self.i];self.i+=1;value=self.atom()
                if c=='_':lo=value
                else:hi=value
            if lo is not None or hi is not None:
                x=el('msubsup' if lo is not None and hi is not None else 'msub' if lo is not None else 'msup')
                x.append(base)
                if lo is not None:x.append(lo)
                if hi is not None:x.append(hi)
                base=x
            a.append(base)
        return a
    def atom(self):
        self.skip()
        if self.i>=len(self.s):raise ValueError('Missing argument: '+self.s)
        c=self.s[self.i];self.i+=1
        if c=='{':return self.seq('}')
        if c=='\\':
            m=re.match(r'[A-Za-z]+',self.s[self.i:])
            if m:
                cmd=m.group();self.i+=len(cmd)
            else:
                cmd=self.s[self.i];self.i+=1
            if cmd in {'left','right','big','Big','bigg'}:return self.atom()
            if cmd in {' ', ',', ';','quad'}:return el('mspace',width='0.20em' if cmd!='quad' else '1em')
            if cmd in {'{','}','|'}:return el('mo',{'|':'‖'}.get(cmd,cmd),stretchy='true')
            if cmd=='frac':
                x=el('mfrac');x.append(self.atom());x.append(self.atom());return x
            if cmd=='sqrt':
                x=el('msqrt');x.append(self.atom());return x
            if cmd in {'widehat','hat','bar','overline','widetilde'}:
                x=el('mover',accent='true');x.append(self.atom());x.append(el('mo',{'widehat':'^','hat':'^','bar':'¯','overline':'¯','widetilde':'~'}[cmd]));return x
            if cmd=='overset':
                over=self.atom();base=self.atom();x=el('mover');x.extend([base,over]);return x
            if cmd in {'mathrm','mathbf','mathcal','mathbb','operatorname'}:
                x=el('mstyle',mathvariant={'mathrm':'normal','mathbf':'bold','mathcal':'script','mathbb':'double-struck','operatorname':'normal'}[cmd]);x.append(self.atom());return x
            if cmd in OPS:return el('mi',cmd,mathvariant='normal')
            if cmd in SYMBOLS:return el('mo' if cmd not in 'alpha beta gamma delta epsilon varepsilon eta theta kappa lambda mu pi rho sigma tau phi varphi chi psi omega Delta Phi Psi Omega Sigma'.split() else 'mi',SYMBOLS[cmd])
            raise ValueError('Unsupported math command '+cmd+' in '+self.s)
        if c.isdigit():return el('mn',c)
        if c.isalpha():return el('mi',c)
        return el('mo',c,stretchy='true' if c in '()[]|' else 'false')

def group_delimiters(node):
    """Use explicit MathML fences so Office preserves unary minus in scripts."""
    for child in list(node):
        group_delimiters(child)
    if E.QName(node).localname != 'mrow':
        return
    children=list(node); out=[]; i=0
    pairs={'(':')','[':']','{':'}','|':'|','‖':'‖'}
    def token(x):
        return x.text if E.QName(x).localname=='mo' else None
    while i<len(children):
        opening=token(children[i])
        if opening not in pairs:
            out.append(children[i]);i+=1;continue
        closing=pairs[opening];depth=1;j=i+1
        while j<len(children):
            value=token(children[j])
            if value==closing:depth-=1
            elif value==opening:depth+=1
            if depth==0:break
            j+=1
        if j==len(children):
            out.append(children[i]);i+=1;continue
        content=el('mrow');content.extend(children[i+1:j]);group_delimiters(content)
        fence=el('mfenced',open=opening,close=closing,separators='')
        fence.append(content);out.append(fence);i=j+1
    node[:]=out

def set_font(run,size=11,bold=None,latin='Times New Roman',east='宋体'):
    run.font.name=latin;run.font.size=Pt(size)
    if bold is not None:run.bold=bold
    rpr=run._element.get_or_add_rPr();rf=rpr.rFonts
    if rf is None:rf=OxmlElement('w:rFonts');rpr.insert(0,rf)
    rf.set(qn('w:eastAsia'),east)
def pformat(p,indent=True,after=5):
    pf=p.paragraph_format;pf.first_line_indent=Pt(22) if indent else Pt(0)
    pf.line_spacing=1.2;pf.space_after=Pt(after);pf.widow_control=True
    if indent:
        p.alignment=WD_ALIGN_PARAGRAPH.JUSTIFY
        pi=p._p.get_or_add_pPr().find(qn('w:ind'))
        if pi is not None:pi.set(qn('w:firstLineChars'),'200')
def cell_margins(cell,top=40,bottom=40,left=80,right=80):
    pr=cell._tc.get_or_add_tcPr();mar=OxmlElement('w:tcMar')
    for k,v in dict(top=top,bottom=bottom,left=left,right=right).items():
        x=OxmlElement('w:'+k);x.set(qn('w:w'),str(v));x.set(qn('w:type'),'dxa');mar.append(x)
    pr.append(mar)
def row_nosplit(row):
    x=OxmlElement('w:cantSplit');row._tr.get_or_add_trPr().append(x)
def borders(table,display=True):
    old=table._tbl.tblPr.find(qn('w:tblBorders'))
    if old is not None:table._tbl.tblPr.remove(old)
    b=OxmlElement('w:tblBorders')
    for n in ('top','left','bottom','right','insideH','insideV'):
        t=OxmlElement('w:'+n);t.set(qn('w:val'),'single' if display and n in ('top','bottom') else 'nil');t.set(qn('w:sz'),'8');t.set(qn('w:color'),'000000');b.append(t)
    table._tbl.tblPr.append(b)
    if display:
        for cell in table.rows[0].cells:
            pr=cell._tc.get_or_add_tcPr();cb=OxmlElement('w:tcBorders')
            line=OxmlElement('w:bottom');line.set(qn('w:val'),'single');line.set(qn('w:sz'),'6');line.set(qn('w:color'),'000000')
            cb.append(line);pr.append(cb)
def build(args):
    import io, zipfile
    here=Path(__file__).resolve().parent
    source=Path(args.source) if args.source else here/'appendix_zh.json'
    output=Path(args.output) if args.output else here/'Appendix_KDRCC_zh.docx'
    template=Path(args.template) if args.template else here/'Appendix_KDRCC_zh.docx'
    if not template.exists(): raise FileNotFoundError('An existing formatted appendix is required as the template.')
    package=template.read_bytes()
    doc=Document(io.BytesIO(package))
    blocks=json.loads(source.read_text(encoding='utf-8'))
    transform=E.XSLT(E.parse(args.xsl))
    def omath(s,size=11):
        root=el('math');root.append(MathParser(s).seq());group_delimiters(root)
        x=copy.deepcopy(transform(root).getroot())
        if E.QName(x).localname=='oMathPara':x=copy.deepcopy(x.find(qn('m:oMath')))
        for r in x.findall('.//'+qn('m:r')):
            pr=r.find(qn('w:rPr'))
            if pr is None:pr=OxmlElement('w:rPr');r.insert(0,pr)
            rf=OxmlElement('w:rFonts');rf.set(qn('w:ascii'),'Cambria Math');rf.set(qn('w:hAnsi'),'Cambria Math');pr.append(rf)
            sz=OxmlElement('w:sz');sz.set(qn('w:val'),str(round(size*2)));pr.append(sz)
        return x
    def findp(prefix):return next(p for p in doc.paragraphs if p.text.startswith(prefix))
    prototypes={
      'p':copy.deepcopy(findp('本附录')._p),
      'h1':copy.deepcopy(findp('S.1 ')._p),
      'h2':copy.deepcopy(findp('S1.1')._p),
      'lead':copy.deepcopy(findp('变量参与度与导数统计量')._p),
      'ref':copy.deepcopy(findp('[1]')._p),
      'caption':copy.deepcopy(findp('表 S1')._p)}
    eqproto=next(copy.deepcopy(t._tbl) for t in doc.tables if len(t.rows)==1 and len(t.columns)==2 and '（S' in t.cell(0,1).text)
    datap=next(t for t in doc.tables if len(t.rows)>1)
    prototypes['headcell']=copy.deepcopy(datap.cell(0,0).paragraphs[0]._p)
    prototypes['cell']=copy.deepcopy(datap.cell(1,0).paragraphs[0]._p)
    body=doc._element.body
    sect=copy.deepcopy(body.sectPr)
    for child in list(body):body.remove(child)
    body.append(sect)
    def clearp(pp):
        for child in list(pp):
            if child.tag!=qn('w:pPr'):pp.remove(child)
    def rich(pp,s,typ='p',size=11):
        proto=prototypes[typ]
        rpr=proto.find('./'+qn('w:r')+'/'+qn('w:rPr'))
        for k,v in enumerate(re.split(r'\$([^$]+)\$',s)):
            if k%2:
                addrun(pp,' ',rpr);pp.append(omath(v,size));addrun(pp,' ',rpr)
            elif v:addrun(pp,v,rpr)
    def addrun(pp,s,rpr):
        r=OxmlElement('w:r')
        if rpr is not None:r.append(copy.deepcopy(rpr))
        t=OxmlElement('w:t');t.set(qn('xml:space'),'preserve');t.text=s;r.append(t);pp.append(r)
    def paragraph(s,typ='p',size=11):
        pp=copy.deepcopy(prototypes[typ]);clearp(pp);rich(pp,s,typ,size);body.insert(len(body)-1,pp);return pp
    def keep(pp):
        pr=pp.find(qn('w:pPr'))
        if pr is None:pr=OxmlElement('w:pPr');pp.insert(0,pr)
        k=pr.find(qn('w:keepNext'))
        if k is None:k=OxmlElement('w:keepNext');pr.append(k)
        k.set(qn('w:val'),'1')
    def spacer():
        pp=OxmlElement('w:p');pr=OxmlElement('w:pPr');pp.append(pr)
        sp=OxmlElement('w:spacing');sp.set(qn('w:after'),'0');sp.set(qn('w:before'),'0');sp.set(qn('w:line'),'20');sp.set(qn('w:lineRule'),'exact');pr.append(sp)
        body.insert(len(body)-1,pp)
    chapter=section='0';counter=0;numbers=[]
    for bi,block in enumerate(blocks):
        typ=block['type']
        if typ in ('h1','h2','lead','p','ref'):
            pp=paragraph(block['text'],typ,10 if typ=='ref' else 11)
            if typ in ('h1','h2','lead') or (bi+1<len(blocks) and blocks[bi+1]['type']=='eq'):keep(pp)
            if typ=='h2':chapter,section=re.match(r'S(\d+)\.(\d+)',block['text']).groups();counter=0
        elif typ=='eq':
            counter+=1;number=f'（S{chapter}-{section}-{counter}）';numbers.append(number)
            table=copy.deepcopy(eqproto);ps=table.findall('./'+qn('w:tr')+'/'+qn('w:tc')+'/'+qn('w:p'))
            for pp in ps:clearp(pp)
            ps[0].append(omath(block['latex'].rstrip().rstrip('.'),10.5))
            r=OxmlElement('w:r');rp=OxmlElement('w:rPr');sz=OxmlElement('w:sz');sz.set(qn('w:val'),'20');rp.append(sz);r.append(rp);t=OxmlElement('w:t');t.text=number;r.append(t);ps[1].append(r)
            pr=ps[1].find(qn('w:pPr'));jc=pr.find(qn('w:jc'))
            if jc is None:jc=OxmlElement('w:jc');pr.append(jc)
            jc.set(qn('w:val'),'right')
            body.insert(len(body)-1,table);spacer()
        elif typ=='table':
            keep(paragraph(block['title'],'caption',10.5))
            vals=[block['headers']]+block['rows'];nc=len(vals[0])
            table=doc.add_table(rows=len(vals),cols=nc);table.alignment=WD_TABLE_ALIGNMENT.CENTER;table.autofit=False
            widths=block.get('widths',[16.8/nc]*nc);widths=[16.8*w/sum(widths) for w in widths]
            for j,w in enumerate(widths):table.columns[j].width=Cm(w)
            for ri,row in enumerate(vals):
                for j,s in enumerate(row):
                    cell=table.cell(ri,j);cell.width=Cm(widths[j]);cell_margins(cell,65,65,85,85);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
                    old=cell.paragraphs[0]._p;pp=copy.deepcopy(prototypes['headcell' if ri==0 else 'cell']);clearp(pp)
                    pr=pp.find(qn('w:pPr'))
                    if pr is not None:
                        for kn in list(pr.findall(qn('w:keepNext'))):pr.remove(kn)
                    rich(pp,s,'headcell' if ri==0 else 'cell',9.5);old.getparent().replace(old,pp)
                    if ri==0:keep(pp)
                row_nosplit(table.rows[ri])
            repeat=OxmlElement('w:tblHeader');table.rows[0]._tr.get_or_add_trPr().append(repeat);borders(table)
            spacer()
    assert len(numbers)==len(set(numbers))
    # Preserve every package part other than the document body, including styles,
    # section settings, headers and footers from the supplied formatted document.
    xml=E.tostring(doc._element,xml_declaration=True,encoding='UTF-8',standalone=True)
    result=io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(package)) as zin,zipfile.ZipFile(result,'w',zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():zout.writestr(info,xml if info.filename=='word/document.xml' else zin.read(info.filename))
    output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(result.getvalue())
    index=Path(args.index) if args.index else here/'equation_index.json'
    index.write_text(json.dumps(numbers,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'output':str(output),'blocks':len(blocks),'display_equations':len(numbers)},ensure_ascii=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--xsl',required=True);ap.add_argument('--output');ap.add_argument('--template');ap.add_argument('--source');ap.add_argument('--index');build(ap.parse_args())
