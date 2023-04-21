import os
import math
import random
import seaborn as sns
from Bio import SeqIO
from ete3 import Tree
from PIL import Image
from itolapi import Itol
import multiprocessing as mp
from PyPDF3.pdf import PageObject
from PyPDF3 import PdfFileWriter, PdfFileReader
from ete3 import TextFace, TreeStyle, NodeStyle


def merge_pdf(pdf_1, pdf_2, margin_size, op_pdf):

    page1 = PdfFileReader(open(pdf_1, "rb"), strict=False).getPage(0)
    page2 = PdfFileReader(open(pdf_2, "rb"), strict=False).getPage(0)

    total_width = page1.mediaBox.upperRight[0] + page2.mediaBox.upperRight[0] + margin_size*3
    total_height = max([page1.mediaBox.upperRight[1], page2.mediaBox.upperRight[1]]) + margin_size*2

    new_page = PageObject.createBlankPage(None, total_width, total_height)
    new_page.mergeTranslatedPage(page1, margin_size, (total_height-margin_size-page1.mediaBox.upperRight[1]))
    new_page.mergeTranslatedPage(page2, (page1.mediaBox.upperRight[0] + margin_size*2), margin_size)

    output = PdfFileWriter()
    output.addPage(new_page)
    output.write(open(op_pdf, "wb"))


merge_pdf('/Users/songweizhi/Desktop/111/l.pdf', '/Users/songweizhi/Desktop/111/r.pdf', 66, '/Users/songweizhi/Desktop/111/combined.pdf')


