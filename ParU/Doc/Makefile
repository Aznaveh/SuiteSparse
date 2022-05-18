
default: paru_user_guide.pdf

include ../../SuiteSparse_config/SuiteSparse_config.mk

paru_user_guide.pdf: paru_user_guide.tex paru_user_guide.bib Makefile
	pdflatex paru_user_guide.tex
	bibtex paru_user_guide
	pdflatex paru_user_guide.tex
	pdflatex paru_user_guide.tex
	- $(RM) -r $(CLEAN) *.out

distclean: purge

clean:
	- $(RM) -r $(CLEAN) 

purge: clean
	- $(RM) -r $(PURGE) *.out paru_user_guide.pdf

