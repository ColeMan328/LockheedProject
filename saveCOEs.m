function save(COEs)
    fid = fopen('COEs.txt', 'w');
    fprintf(fid, 'a = %s\n', COEs{1});
    fprintf(fid, 'e = %s\n', COEs{2});
    fprintf(fid, 'i = %s\n', COEs{3});
    fprintf(fid, '\u03C9. = %s\n', COEs{4});
    fprintf(fid, '\u03A9. = %s\n', COEs{5});
    fprintf(fid, '\u03BD. = %s\n', COEs{6});
    fprintf(fid, 'u = %s\n', COEs{7});
    fprintf(fid, '\u03A0. = %s\n', COEs{8});
    fprintf(fid, 'l = %s\n', COEs{9});
    fprintf(fid, 'Orbit = %s\n', COEs{10});

end