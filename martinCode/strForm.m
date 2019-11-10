function stringEq=strForm(eq)

    stringEq=sprintf('%s',eq);
    stringEq=regexprep(stringEq,'Kx','kx');
    stringEq=regexprep(stringEq,'c(\d)1','c$1(pl)');
    stringEq=regexprep(stringEq,'c(\d)2','c$1(il)');
    stringEq=regexprep(stringEq,'c(\d)3','c$1(ol)');
    stringEq=regexprep(stringEq,'kx1','kx(pl)');
    stringEq=regexprep(stringEq,'kx2','kx(il)');
    stringEq=regexprep(stringEq,'kx3','kx(ol)');
    stringEq=regexprep(stringEq,'*','.*');
    stringEq=regexprep(stringEq,'/','./');

end