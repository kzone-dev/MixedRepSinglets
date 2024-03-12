using Pkg; Pkg.activate(".")
using DelimitedFiles
using LaTeXStrings

function write_tex_table(name,data;insert_hline=[],no_header=false)
    io = open(name,"w")
    rows, cols = size(data)
    # I have hardcoded centering of the table's contents
    table_layout = repeat("|c",cols)*"|"
    header = """\\begin{tabular}{$table_layout}
    \t\\hline
    """
    # either write a header and start with the second row or print the header
    write(io,header)
    for i in 1:rows
        # if the gauge coupling changes insert two hlines for formatting
        if i ∈ insert_hline
            write(io,"\t\\hline\\hline\n")
        end
        # we could use some padding here for nicer formatting but for now
        # I only insert a '&' to make a minimal tex-compliant table
        write(io,"\t"*string(data[i,1])*"&")
        for j in 2:cols-1
            write(io,string(data[i,j])*"&")
        end
        write(io,string(data[i,cols])*"\\\\\n")
        !no_header && i==1 && write(io,"\t\\hline\\hline\n")
    end
    write(io,"\t\\hline\\hline\n")
    write(io,"\\end{tabular}")
    close(io)
end

results = readdlm("output/tables/table_results.csv",';')
fitting = readdlm("output/tables/table_fitting.csv",';')
gevp = readdlm("output/tables/table_gevp.csv",';')

write_tex_table("output/tables/table_results.tex",results)
write_tex_table("output/tables/table_fitting.tex",fitting)
write_tex_table("output/tables/table_gevp.tex",gevp)
