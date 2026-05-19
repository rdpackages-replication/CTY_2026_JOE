/**************************************************************************************************
* Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods
* Empirical Application: SPP Data
* Stata replication output
**************************************************************************************************/

clear all
set more off
set seed 20260501

capture confirm file "spp.csv"
if _rc {
    local dofile = subinstr("`c(filename)'", "\\", "/", .)
    local slash = strrpos("`dofile'", "/")
    if `slash' > 0 cd "`=substr("`dofile'", 1, `slash' - 1)'"
}

local rd2d_stata_path : environment RD2D_STATA_PATH
if "`rd2d_stata_path'" != "" adopath ++ "`rd2d_stata_path'"
capture which rd2d_dist
if _rc {
    di as error "Stata rd2d package not found. Install rd2d or set RD2D_STATA_PATH to the rd2d/stata folder."
    exit 111
}

local outdir : environment RD2D_OUTPUT_DIR
if "`outdir'" == "" local outdir "output"
cap mkdir "`outdir'"
local oldfiles : dir "`outdir'" files "empapp_*.csv"
foreach f of local oldfiles {
    erase "`outdir'/`f'"
}

local refdir : environment RD2D_REFERENCE_OUTPUT
local repp : environment RD2D_EMP_REPP
if "`repp'" == "" local repp 5000

mata:
string scalar rd2d_emp_fmt(real scalar x)
{
    if (x >= .) return("")
    return(strofreal(x, "%21.15g"))
}

real scalar rd2d_emp_mean(real colvector x)
{
    x = select(x, x :< .)
    if (rows(x) == 0) return(.)
    return(mean(x))
}

real scalar rd2d_emp_max(real colvector x)
{
    x = select(x, x :< .)
    if (rows(x) == 0) return(.)
    return(max(x))
}

string scalar rd2d_emp_row(string scalar label, real rowvector values)
{
    real scalar j
    string scalar line
    line = label
    for (j=1; j<=cols(values); j++) line = line + "," + rd2d_emp_fmt(values[j])
    return(line)
}

void rd2d_emp_add_wbate(real matrix V, real scalar level, real rowvector values)
{
    real scalar n, se, se2, center, tval, cval
    real rowvector w
    n = rows(V)
    if (n > 0 & rows(V) == cols(V)) {
        w = J(1, n, 1/n)
        se2 = (w * V * w')[1,1]
        if (se2 < 0) se2 = 0
        se = sqrt(se2)
        center = values[12]
        if (se > 0 & center < .) {
            tval = center / se
            cval = invnormal((level + 100) / 200)
            values[16] = se
            values[13] = tval
            values[14] = 2 * (1 - normal(abs(tval)))
            values[8] = center - cval * se
            values[9] = center + cval * se
        }
    }
}

void rd2d_emp_write_metadata(string scalar file, string scalar repp)
{
    real scalar fh
    fh = fopen(file, "w")
    fput(fh, "name,value")
    fput(fh, "rd2d.stata,0.1.0")
    fput(fh, "repp," + repp)
    fput(fh, "methods,smooth; adaptive; unknown_kink; rdrobust")
    fclose(fh)
}

void rd2d_emp_write_distance_csv(string scalar file, string scalar mname, string scalar vname)
{
    real matrix M, V
    real scalar fh, i
    real rowvector vals, wb, lb
    M = st_matrix(mname)
    V = (vname == "" ? J(0,0,.) : st_matrix(vname))
    fh = fopen(file, "w")
    fput(fh, "row,b1,b2,h0,h1,N.Co,N.Tr,estimate.p,ci.lower,ci.upper,cb.lower,cb.upper,estimate.q,t.value,p.value,std.err.p,std.err.q,h0.rbc,h1.rbc")
    for (i=1; i<=rows(M); i++) {
        vals = (M[i,1], M[i,2], M[i,11], M[i,12], M[i,15], M[i,16], M[i,3], M[i,9], M[i,10], ., ., M[i,5], M[i,7], M[i,8], M[i,4], M[i,6], M[i,13], M[i,14])
        fput(fh, rd2d_emp_row(strofreal(i), vals))
    }
    wb = J(1, 18, .)
    wb[7] = rd2d_emp_mean(M[,3])
    wb[12] = rd2d_emp_mean(M[,5])
    rd2d_emp_add_wbate(V, 95, wb)
    fput(fh, rd2d_emp_row("WBATE", wb))
    lb = J(1, 18, .)
    lb[7] = rd2d_emp_max(M[,3])
    lb[12] = rd2d_emp_max(M[,5])
    fput(fh, rd2d_emp_row("LBATE", lb))
    fclose(fh)
}
end

program define rd2d_emp_numlist_from_csv, rclass
    syntax using/, Columns(string)
    capture confirm file "`using'"
    if _rc {
        return local list ""
        exit
    }
    preserve
    import delimited using "`using'", clear varnames(1)
    capture confirm variable row
    if !_rc {
        gen double __rownum = real(row)
        keep if __rownum < .
    }
    local values ""
    forvalues i = 1/`=_N' {
        foreach col of local columns {
            local values `"`values' `=string(`col'[`i'], "%21.15g")'"'
        }
    }
    restore
    return local list `"`values'"'
end

import delimited using "spp.csv", clear varnames(1)
keep running_saber11 running_sisben spadies_any eligible_spp beneficiary_spp
rename running_saber11 x1
rename running_sisben x2
rename spadies_any y
rename eligible_spp d
rename beneficiary_spp w
drop if missing(x1, x2, y, d, w)
gen byte __expected_assignment = x1 >= 0 & x2 >= 0
assert d == __expected_assignment
drop __expected_assignment
frame rename default spp

frame create eval
frame eval {
    set obs 40
    gen double x1 = cond(_n <= 20, 0, (_n - 21) * 56 / 20)
    gen double x2 = cond(_n <= 20, 40 - (_n - 1) * 40 / 20, 0)
    keep in 11/31
}

frame spp {
    quietly summarize x1
    local sx1 = r(sd)
    quietly summarize x2
    local sx2 = r(sd)
    local scale2 = `sx1' / `sx2'
    replace x2 = x2 * `scale2'
}
frame eval {
    replace x2 = x2 * `scale2'
    local blist ""
    forvalues i = 1/`=_N' {
        local bx`i' = x1[`i']
        local by`i' = x2[`i']
        local blist `"`blist' `=string(x1[`i'], "%21.15g")' `=string(x2[`i'], "%21.15g")'"'
    }
}

frame change spp
local distvars ""
forvalues j = 1/21 {
    gen double dist`j' = sqrt((x1 - `bx`j'')^2 + (x2 - `by`j'')^2) * (2*d - 1)
    local distvars "`distvars' dist`j'"
}

mata: rd2d_emp_write_metadata("`outdir'/empapp_metadata.csv", "`repp'")

foreach method in smooth adaptive unknown_kink rdrobust {
    local hopt ""
    if "`refdir'" != "" {
        rd2d_emp_numlist_from_csv using "`refdir'/empapp_`method'_itt.csv", columns(h0 h1)
        if "`r(list)'" != "" local hopt "h(`r(list)')"
    }
    if "`method'" == "rdrobust" & "`hopt'" == "" {
        capture which rdbwselect
        if _rc {
            di as error "rdrobust/rdbwselect is needed for the rdrobust bandwidth rule unless RD2D_REFERENCE_OUTPUT is set."
            exit 111
        }
        local hlist ""
        forvalues j = 1/21 {
            quietly rdbwselect y dist`j', vce(hc1)
            matrix __bws = e(mat_h)
            local hlist `"`hlist' `=string(__bws[1,1], "%21.15g")' `=string(__bws[1,2], "%21.15g")'"'
        }
        local hopt "h(`hlist')"
    }

    local kinkopt ""
    if "`method'" == "unknown_kink" {
        if "`hopt'" == "" local kinkopt "kinkunknown(1 0)"
        else local kinkopt "q(1)"
    }
    if "`method'" == "adaptive" & "`hopt'" == "" local kinkopt "kinkposition(11)"

    quietly rd2d_dist y `distvars', b(`blist') vce(hc1) repp(`repp') `kinkopt' `hopt'
    tempname sharp
    matrix `sharp' = e(main)

    quietly rd2d_dist y `distvars', b(`blist') fuzzy(w) bwparam(itt) paramscov(main itt fs) vce(hc1) repp(`repp') `kinkopt' `hopt'
    matrix __diff = `sharp' - e(itt)
    mata: st_numscalar("__maxdiff", max(abs(st_matrix("__diff"))))
    if scalar(__maxdiff) > 1e-8 {
        di as error "Sharp main and fuzzy ITT differ for `method'."
        exit 498
    }
    mata: rd2d_emp_write_distance_csv("`outdir'/empapp_`method'_fuzzy.csv", "e(main)", "e(V_main)")
    mata: rd2d_emp_write_distance_csv("`outdir'/empapp_`method'_itt.csv", "e(itt)", "e(V_itt)")
    mata: rd2d_emp_write_distance_csv("`outdir'/empapp_`method'_fs.csv", "e(fs)", "e(V_fs)")
}

di as text "Stata empirical application output data complete."
di as text "CSV outputs written to: `outdir'"


