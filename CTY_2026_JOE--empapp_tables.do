/**************************************************************************************************
* Distance-Based Methods: empirical application tables from Stata CSV output
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
local tablesdir : environment RD2D_TABLES_DIR
if "`tablesdir'" == "" local tablesdir "tables"
cap mkdir "`tablesdir'"
local oldfiles : dir "`tablesdir'" files "empapp_*.tex"
foreach f of local oldfiles {
    erase "`tablesdir'/`f'"
}

program define write_distance_table
    syntax using/, Saving(string)
    preserve
    import delimited using "`using'", clear varnames(1)
    gen double __rownum = real(row)
    keep if inlist(row,"WBATE","LBATE") | __rownum < .
    file open tab using "`saving'", write replace text
    file write tab "\begin{tabular}{lrrrrrc}" _n
    file write tab "\toprule" _n
    file write tab "Boundary & h & N.Co & N.Tr & Estimate & p-value & 95\% CI \\\" _n
    file write tab "\midrule" _n
    forvalues i=1/`=_N' {
        local label = row[`i']
        if inlist("`label'", "WBATE", "LBATE") file write tab "\midrule" _n
        local h = cond(missing(h0[`i']), "", string(h0[`i'], "%9.3f"))
        local n0 = cond(missing(nco[`i']), "", string(nco[`i'], "%9.0f"))
        local n1 = cond(missing(ntr[`i']), "", string(ntr[`i'], "%9.0f"))
        local est = cond(missing(estimatep[`i']), "", string(estimatep[`i'], "%9.3f"))
        local pv = cond(missing(pvalue[`i']), "", string(pvalue[`i'], "%9.3f"))
        local ci = cond(missing(cilower[`i']) | missing(ciupper[`i']), "", "(" + string(cilower[`i'], "%9.3f") + ", " + string(ciupper[`i'], "%9.3f") + ")")
        file write tab "`label' & `h' & `n0' & `n1' & `est' & `pv' & `ci' \\\" _n
    }
    file write tab "\bottomrule" _n
    file write tab "\end{tabular}" _n
    file close tab
    restore
end

foreach method in smooth adaptive unknown_kink rdrobust {
    foreach estimand in fuzzy itt fs {
        write_distance_table using "`outdir'/empapp_`method'_`estimand'.csv", saving("`tablesdir'/empapp_`method'_`estimand'.tex")
    }
}

di as text "Stata empirical application tables written to: `tablesdir'"


