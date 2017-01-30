latexfile = IKthesis
TEX = latex -halt-on-error
#BIB = latex -halt-on-error
BIB= bibtex
MPOST = mpost
release: $(latexfile).pdf
	cp $(latexfile).pdf $(latexfile)_svn$(shell svn info | grep Revision| cut -f 2 -d\ ).pdf

all: $(latexfile).pdf
#IKthesis.bst: IKthesis.dbj
#	$(TEX) IKthesis.dbj


$(latexfile).dvi : IKthesis.bst IKthesis-ack-txt.tex IKthesis-stm-tab.tex IKthesis-stm-fig.tex    IKthesis-stm-txt.tex  IKthesis.tex IKthesis-tit.tex
	sed -i  's/\draftversion.*/\draftversion{svn'$$(svn info | grep Revision| cut -f 2 -d\ )'}/g'  IKthesis-tit.tex
	$(TEX) $(latexfile)
	for a in $$(cd Metapost; ls -1 *mp);\
	do \
	cd Metapost;\
	$(MPOST) $$a;\
	done
	$(BIB) $(latexfile)
	$(BIB)  $(latexfile)
	$(TEX) $(latexfile)
	$(TEX) $(latexfile)

$(latexfile).pdf : $(latexfile).ps
	ps2pdf -dPDFSETTINGS=/printer -sPAPERSIZE=a4  $(latexfile).ps $(latexfile).pdf

pdf : $(latexfile).pdf
$(latexfile).ps : $(latexfile).dvi
	dvips -Ppdf -t a4 $(latexfile)

ps : $(latexfile).ps 




allpres: Presentation/presentation_$(latexfile).pdf

presentation_$(latexfile).dvi :  Presentation/presentation_$(latexfile).tex
	$(TEX) Presentation/presentation_$(latexfile)
#	mv presentation_$(latexfile)* Presentation/
#	$(TEX) $(latexfile)
	for a in $$(cd Metapost; ls -1 *mp);\
	do \
	cd Metapost;\
	$(MPOST) $$a;\
	done
#	bibtex Presentation/presentation_$(latexfile) 
#	cp $(latexfile).bbl Presentation/presentation_$(latexfile).bbl 
#	cp $(latexfile).bst Presentation/presentation_$(latexfile).bst
#	mv presentation_$(latexfile)* Presentation/
	$(TEX) Presentation/presentation_$(latexfile)
#	mv presentation_$(latexfile)* Presentation/
	$(TEX) Presentation/presentation_$(latexfile)
#	mv presentation_$(latexfile)* Presentation/


Presentation/presentation_$(latexfile).pdf : Presentation/presentation_$(latexfile).ps
	ps2pdf -dPDFSETTINGS=/printer -sPAPERSIZE=a4 Presentation/presentation_$(latexfile).ps Presentation/presentation_$(latexfile).pdf

pdfpres : Presentation/presentation_$(latexfile).pdf

Presentation/presentation_$(latexfile).ps : presentation_$(latexfile).dvi
	dvips  presentation_$(latexfile).dvi
	mv  presentation_$(latexfile).* Presentation

pspres : Presentation/presentation_$(latexfile).ps 





clean :
	rm -f Metapost/*
	rm -f *.aux
	rm -f *.lof
	rm -f *.out
	rm -f *.lot
	rm -f *.bbl
#	rm -f *.bst
	rm -f *.l4l
	rm -f *.log
	rm -f *.blg
	rm -f *.ilg
	rm -f *.aux
	rm -f *.dvi
	rm -f *.ps
	rm -f *.pdf
	rm -f *.toc
	rm -f *.idx
	rm -f *.ind
	rm -f *.nav
	rm -f *.snm
	rm -f Presentation/*.aux
	rm -f Presentation/*.lof
	rm -f Presentation/*.out
	rm -f Presentation/*.lot
	rm -f Presentation/*.bbl
	rm -f Presentation/*.l4l
	rm -f Presentation/*.log
	rm -f Presentation/*.blg
	rm -f Presentation/*.ilg
	rm -f Presentation/*.aux
	rm -f Presentation/*.dvi
	rm -f Presentation/*.ps
	rm -f Presentation/*.pdf
	rm -f Presentation/*.toc
	rm -f Presentation/*.idx
	rm -f Presentation/*.ind
	rm -f Presentation/*.nav
	rm -f Presentation/*.snm

