@default_files = ("tcc.tex");
$pdf_mode = 1;
$dvi_previewer = "start evince %O %S";
$pdf_previewer = "start evince %O %S";
$dvi_update_method=0;
$pdf_update_method=0;
$pdflatex = 'pdflatex  %O  --shell-escape %S';

# Custom dependency and function for nomencl package 
add_cus_dep( 'nlo', 'nls', 0, 'makenlo2nls' );
sub makenlo2nls {
    system( "makeindex -s nomencl.ist -o \"$_[0].nls\" \"$_[0].nlo\"" );
}
