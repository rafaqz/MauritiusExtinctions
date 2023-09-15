for pdf in readdir("/home/raf/PhD/Mascarenes/maps/Mauritius/Studies_of_Mascarine_birds"; join=true)
    name, ext = splitext(pdf)
    ext == ".pdf" || continue
    png = name * ".png" 
    run(`pdftoppm $pdf $png -png -r 300`)
end
