# makefile for analysis pdfs
PDFS=

pdf: $(PDFS)

%.pdf: %.inpynb
	jupyter nbconvert $< --to pdf
	