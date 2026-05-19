/**************************************************************************************************
* Distance-Based Methods: empirical application figures from Stata CSV output
**************************************************************************************************/
clear all
set more off

capture confirm file "spp.csv"
if _rc {
    local dofile = subinstr("`c(filename)'", "\\", "/", .)
    local slash = strrpos("`dofile'", "/")
    if `slash' > 0 cd "`=substr("`dofile'", 1, `slash' - 1)'"
}

local outdir : environment RD2D_OUTPUT_DIR
if "`outdir'" == "" local outdir "output"
local figuresdir : environment RD2D_FIGURES_DIR
if "`figuresdir'" == "" local figuresdir "figures"
cap mkdir "`figuresdir'"

program define effect_plot
    syntax using/, Saving(string) Ytitle(string)
    preserve
    import delimited using "`using'", clear varnames(1)
    gen double idx = real(row)
    keep if idx < .
    twoway (rcap cilower ciupper idx, lcolor(gs10)) ///
           (scatter estimatep idx, mcolor(black) msize(vsmall)), ///
           legend(off) xtitle("Boundary point") ytitle("`ytitle'") graphregion(color(white))
    graph export "`saving'", replace width(1800)
    restore
end

foreach method in smooth adaptive unknown_kink rdrobust {
    effect_plot using "`outdir'/empapp_`method'_itt.csv", saving("`figuresdir'/fig4-`method'.png") ytitle("Treatment effect")
}

di as text "Stata empirical application figures written to: `figuresdir'"

