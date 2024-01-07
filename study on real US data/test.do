import excel "/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne/2021data_stata.xlsx", sheet("Sheet1") firstrow clear
gen times=quarterly(date,"YQ")	
tsset times, quarterly

clear matrix
set matsize 1000

cap drop y p c i h
gen y=log(RealGDP)
gen p=log(CPI)
gen consumption = RealConsumptionNonDurables + RealConsumptionService
gen c = log(consumption)
gen investment = RealConsumptionDurables + RealInvestment
gen i = log(investment)
gen h = log(TotalHours)
rename ShadowRate s
rename TFP t
gen trend=_n

drop if time<1980

******************
*  Unit Root tests: only i and t are stationay
******************

// scalar T = 247
scalar pmax = int(12*((247+1)/100)^0.25)
display pmax

/* testing y */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg y trend l.y l(1/15)d.y
qui estat ic
matrix yorder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg y trend l.y l(1/`i')d.y if _n>preg_l
qui estat ic
matrix yorder_l=(yorder_l\r(S))
}

qui reg y trend l.y if _n>preg_l
qui estat ic
matrix yorder_l=(yorder_l\r(S))

matlist yorder_l     /* AIC: P=3, SIC: P=3 */
dfuller y, trend lag(0)  

/* testing d.y */
cap drop diff
scalar diff = 1
scalar preg_d = pmax+1+diff
scalar list pmax diff preg_d

qui reg d.y ld.y l(1/15)d2.y
qui estat ic
matrix yorder_d=r(S)

forval j=1/15 {
local k=15-`j'
qui reg d.y ld.y l(1/`k')d2.y if _n>preg_d
qui estat ic
matrix yorder_d=(yorder_d\r(S))
}

qui reg d.y ld.y if _n>preg_d 
qui estat ic
matrix yorder_d=(yorder_d\r(S))

matlist yorder_d     /* AIC: P=3, SIC: P=3 */
dfuller d.y, lag(1)  /* I(1) */

/* testing p */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg p trend l.p l(1/15)d.p
qui estat ic
matrix porder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg p trend l.p l(1/`i')d.p if _n>preg_l
qui estat ic
matrix porder_l=(porder_l\r(S))
}

qui reg p trend l.p if _n>preg_l
qui estat ic
matrix porder_l=(porder_l\r(S))

matlist porder_l     /* AIC: P=5, SIC: P=3 */
dfuller p, trend lag(6)  
dfuller p, trend lag(3)  

/* testing d.p */
cap drop diff
scalar diff = 1
scalar preg_d = pmax+1+diff
scalar list pmax diff preg_d

qui reg d.p ld.p l(1/15)d2.p
qui estat ic
matrix porder_d=r(S)

forval j=1/15 {
local k=15-`j'
qui reg d.p ld.p l(1/`k')d2.p if _n>preg_d
qui estat ic
matrix porder_d=(porder_d\r(S))
}

qui reg d.p ld.p if _n>preg_d 
qui estat ic
matrix porder_d=(porder_d\r(S))

matlist porder_d     /* AIC: P=4, SIC: P=2 */
dfuller d.p, lag(4)  /* I(1) */
dfuller d.p, lag(1)  /* I(1) */

/* testing c */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg c trend l.c l(1/15)d.c
qui estat ic
matrix corder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg c trend l.c l(1/`i')d.c if _n>preg_l
qui estat ic
matrix corder_l=(corder_l\r(S))
}

qui reg c trend l.c if _n>preg_l
qui estat ic
matrix corder_l=(corder_l\r(S))

matlist corder_l     /* AIC: P=1, SIC: P=4 */
dfuller c, trend lag(1)  
dfuller c, trend lag(4)  

/* testing d.c */
cap drop diff
scalar diff = 1
scalar preg_d = pmax+1+diff
scalar list pmax diff preg_d

qui reg d.c ld.c l(1/15)d2.c
qui estat ic
matrix corder_d=r(S)

forval j=1/15 {
local k=15-`j'
qui reg d.c ld.c l(1/`k')d2.c if _n>preg_d
qui estat ic
matrix corder_d=(corder_d\r(S))
}

qui reg d.c ld.c if _n>preg_d 
qui estat ic
matrix corder_d=(corder_d\r(S))

matlist corder_d     /* AIC: P=0, SIC: P=0 */
dfuller d.c, lag(0)  /* I(1) */

/* testing i */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg i trend l.i l(1/15)d.i
qui estat ic
matrix iorder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg i trend l.i l(1/`i')d.i if _n>preg_l
qui estat ic
matrix iorder_l=(iorder_l\r(S))
}

qui reg i trend l.i if _n>preg_l
qui estat ic
matrix iorder_l=(iorder_l\r(S))

matlist iorder_l     /* AIC: P=2, SIC: P=1 */
dfuller i, trend lag(2)  
dfuller i, trend lag(1)  /* stationary */

/* testing h */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg h trend l.h l(1/15)d.h
qui estat ic
matrix horder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg h trend l.h l(1/`i')d.h if _n>preg_l
qui estat ic
matrix horder_l=(horder_l\r(S))
}

qui reg h trend l.h if _n>preg_l
qui estat ic
matrix horder_l=(horder_l\r(S))

matlist horder_l     /* AIC: P=0, SIC: P=0 */
dfuller h, trend lag(0)  

/* testing d.h */
cap drop diff
scalar diff = 1
scalar preg_d = pmax+1+diff
scalar list pmax diff preg_d

qui reg d.h ld.h l(1/15)d2.h
qui estat ic
matrix horder_d=r(S)

forval j=1/15 {
local k=15-`j'
qui reg d.h ld.h l(1/`k')d2.h if _n>preg_d
qui estat ic
matrix horder_d=(horder_d\r(S))
}

qui reg d.h ld.h if _n>preg_d 
qui estat ic
matrix horder_d=(horder_d\r(S))

matlist horder_d     /* AIC: P=0, SIC: P=0 */
dfuller d.h, lag(0)  /* I(1) */

/* testing s */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg s trend l.s l(1/15)d.s
qui estat ic
matrix sorder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg s trend l.s l(1/`i')d.s if _n>preg_l
qui estat ic
matrix sorder_l=(sorder_l\r(S))
}

qui reg s trend l.s if _n>preg_l
qui estat ic
matrix sorder_l=(sorder_l\r(S))

matlist sorder_l     /* AIC: P=7, SIC: P=3 */
dfuller s, trend lag(7)  
dfuller s, trend lag(3)  

/* testing d.s */
cap drop diff
scalar diff = 1
scalar preg_d = pmax+1+diff
scalar list pmax diff preg_d

qui reg d.s ld.s l(1/15)d2.s
qui estat ic
matrix sorder_d=r(S)

forval j=1/15 {
local k=15-`j'
qui reg d.s ld.s l(1/`k')d2.s if _n>preg_d
qui estat ic
matrix sorder_d=(sorder_d\r(S))
}

qui reg d.s ld.s if _n>preg_d 
qui estat ic
matrix sorder_d=(sorder_d\r(S))

matlist sorder_d     /* AIC: P=6, SIC: P=2 */
dfuller d.s, lag(6)  /* I(1) */
dfuller d.s, lag(2)  /* I(1) */

/* testing t */
scalar diff = 0
scalar preg_l = pmax + 1 + diff
scalar list pmax diff preg_l 

qui reg t trend l.t l(1/15)d.t
qui estat ic
matrix torder_l=r(S)

forval j=1/14 {
local i=15-`j'
qui reg t trend l.t l(1/`i')d.t if _n>preg_l
qui estat ic
matrix torder_l=(torder_l\r(S))
}

qui reg t trend l.t if _n>preg_l
qui estat ic
matrix torder_l=(torder_l\r(S))

matlist torder_l     /* AIC: P=0, SIC: P=0 */
dfuller t, trend lag(0)  /* stationary */

********
*  VAR *
********
varsoc d.y d.p d.c i d.h d.s t, maxlag(8)   /* p = 1*/
varsoc y p i c h s t, maxlag(8)   /* p = 1, 2 */

/* estimation and diagnosis_level case */
var y p i c h s t, lag(1)
varlmar, mlag(8)     /* serially correlated */
varstable            /* roots are fine */
varnorm              /* D_y residual not normal, kurtosis issue */

irf set oir_ovd, replace
qui var y p i c h s t, lag(1/8)
irf create order1, step(51) replace
irf graph oirf, irf(order1) impulse(y p) response(y p) level(95)

/* estimation and diagnosis_diff case */
var d.y d.p d.c i d.h d.s t, lag(1)
varlmar, mlag(8)     /* serially correlated */
varstable            /* roots are fine */
varnorm  

gen dy = d.y
gen dp = d.p
irf set oir_ovd_, replace
qui var dy dp d.c i d.h d.s t, lag(1)
irf create order1, step(51) replace
irf graph oirf, irf(order1) impulse(dy dp) response(dy dp) level(95)
