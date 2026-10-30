"""Post-process pandoc docx: figure width, centred figures, table borders, column widths, 9-pt table text."""
import docx, sys, zipfile, shutil, re
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


def fix(p, out, total=6.7):
    d = docx.Document(p)
    for sh in d.inline_shapes:
        if sh.width > Inches(6.5) or sh.width < Inches(3.0):
            r = sh.height / sh.width; w = Inches(6.5) if sh.width > Inches(6.5) else sh.width; sh.width = w; sh.height = int(w * r)
    for para in d.paragraphs:
        if para._p.xpath('.//w:drawing'):
            para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for t in d.tables:
        tblPr = t._tbl.tblPr
        for old in tblPr.findall(qn('w:tblBorders')):
            tblPr.remove(old)
        b = OxmlElement('w:tblBorders')
        for e in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
            x = OxmlElement('w:' + e); x.set(qn('w:val'), 'single'); x.set(qn('w:sz'), '4'); x.set(qn('w:space'), '0'); x.set(qn('w:color'), '808080'); b.append(x)
        anchor = None
        for tag in ['tblInd', 'tblCellSpacing', 'jc', 'tblW', 'tblStyleColBandSize', 'tblStyleRowBandSize', 'bidiVisual', 'tblOverlap', 'tblpPr', 'tblStyle']:
            el = tblPr.find(qn('w:' + tag))
            if el is not None:
                anchor = el; break
        anchor.addnext(b) if anchor is not None else tblPr.insert(0, b)
        n = len(t.columns)
        first = 1.7 if n > 2 else total / 2
        widths = [first] + [(total - first) / (n - 1)] * (n - 1)
        for row in t.rows:
            for i, c in enumerate(row.cells):
                c.width = Inches(widths[i])
                for pp in c.paragraphs:
                    for rr in pp.runs:
                        rr.font.size = Pt(9)
        for c in t.rows[0].cells:
            for pp in c.paragraphs:
                for rr in pp.runs:
                    rr.bold = True
    d.core_properties.author = ''
    d.save(out)
    zin = zipfile.ZipFile(out); tmp = out + '.tmp'; zo = zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED)
    for it in zin.infolist():
        bb = zin.read(it.filename)
        if it.filename == 'word/document.xml':
            bb = re.sub(rb'w:w="(\d+)\.0"', rb'w:w="\1"', bb)
        zo.writestr(it, bb)
    zo.close(); zin.close(); shutil.move(tmp, out)


if __name__ == '__main__':
    fix(sys.argv[1], sys.argv[2])
