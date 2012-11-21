rm -f mg_matrix.o mg_matrix.so
R CMD SHLIB mg_matrix.c -lRlapack
strip -s mg_matrix.so
