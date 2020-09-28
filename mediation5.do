**************************************************************************
*
* program to simulate continuous outcome data (change in pain intensity) for 5 parallel mediators
* for grant with Mark Jensen and Melissa Day
*
* Author: Susanne May
* Date: 09/19/2020
*
*
* this program requires input parameters, which need to be supplied
*
* These input parameters are: 
*
*        1 = nperg         number of observations per group for two groups (placebo and MM for now)
*        2 = toteff        total (direct and indirect) effect of the intervention on change in pain intensity (assumed standardized to SD=1)
*        3 = indpoftot     indirect is percent of total (e.g. 0.5 for 50% or half of the total effect represents the indirect effect)
*        4 = toteffsd      SD of total effect (assume 1 for now)
*        5 = seed        number seed
*        6 = reps          number of replications 
*
*       . do mediation5 50 -.40 .5 1 979083 100
*
* run this file once as 
*  
*    do mediation5
*
* before running runlots.do
*
**************************************************************************
/*
*** for singlerun
  local nperg=50

  local toteff=-0.40
  local indpoftot=0.50
  local toteffsd=1
  
    
  local seed=528320
  local reps=1
*  local no=01
  
  set more off
  set seed `seed'

*** end for singlerun
*/
clear

*capture log close
*** for single run use next line
*capture log using mediation5-onerun.log, replace

*** comment out next few lines for single run
capture program drop simula
program define simula, rclass

        version 16.0
        
        *----------------
        syntax [, nperg(integer 1) toteff(real 1) indpoftot(real 1) toteffsd(real 1) seed(integer 1) ]
        *----------------

        drop _all
        
        display "nperg="`nperg' "  toteff=" `toteff' "  indpoftot=="`indpoftot' "  toteffsd=" `toteffsd' "  seed=" `seed'
        
        * create counter for potentially looking at individual simulated data sets
        use one
        local nummerin=nummer[1]
        local nummerout=`nummerin'+1
        replace nummer=`nummerout' in 1
        save,replace
        drop _all

        local total=2*`nperg'
        
        set obs `total'
        gen tx=(_n>=`total'/2)
        
        local indeff=`toteff'*`indpoftot'
        local direff=`toteff'-`indeff'
        display "toteff = " `toteff' "  indpoftot = " `indpoftot' "  indeff = " `indeff' "  direff = " `direff' 
        local indab=`indeff'/5
        
        *** for now assume that each of the factors has the same indirect effect in magnitude and that a and b are the same in magnitude (thus the sqrt)
        local ai=sqrt(abs(`indab'))
        local bi=sign(`indab')*sqrt(abs(`indab'))
        display "   ai = " `ai' "   bi = " `bi' "   indab = " `indab'
        
        *** modeling the mediators on the tx 
        gen M1=invnorm(uniform())+`ai'*tx
        gen M2=invnorm(uniform())+`ai'*tx
        gen M3=invnorm(uniform())+`ai'*tx
        gen M4=invnorm(uniform())+`ai'*tx
        gen M5=invnorm(uniform())+`ai'*tx
        
        gen deltaY=invnorm(uniform())+`bi'*M1+`bi'*M2+`bi'*M3+`bi'*M4+`bi'*M5+`direff'*tx
        
        sum


 /*
 **** next section is left over from previous program, but left in if I decide to model the mediators with correlation 
 
        * create appropriate bivariate normally distributed outcomes for the WL group (no effect) 
        * create a vector that contains the equivalent of a lower triangular correlation matrix
        matrix cthetawl = (1, `corrtheta', 1, 0, 0, 1, 0, 0, `corrfatig', 1)
*        matrix cthetawl = (1, `corrtheta', 0, 0 \ `corrtheta', 1, 0, 0 \ 0, 0, 1, `corrfatig' \ 0, 0, `corrfatig', 1)

        * create a vector that contains the means of the variables for the waitlist group
	    matrix mthetawl = (`meanthetabl',`meanthetabl',`meanfatigbl', `meanfatigbl')

	    * create a vector that contains the standard deviations
	    matrix sdthetawl = (`sdthetabl',`sdthetabl',`sdfatigbl',`sdfatigbl')
	    * draw a sample of specified cases from a normal distribution with specified correlation structure
	    * and specified means and standard deviations

	    drawnorm thetabl thetapost fatigbl fatigpost, n(`nperg') corr(cthetawl) cstorage(lower) means(mthetawl) sds(sdthetawl) double
        gen tx=0
        save datawl, replace
        drop _all

        * create appropriate bivariate normally distributed outcomes for the NF+HYP group (effect under alternative) 
        * create a vector that contains the equivalent of a lower triangular correlation matrix
        matrix cthetanfhyp = (1, `corrtheta', 1, 0, 0, 1, 0, 0, `corrfatig', 1)
        
	    * create a vector that contains the means of the variables for the treatment group
        ************
        ************ this is the main code for generating the treatment effect (for now without the mediation component)
        ************
	    matrix mthetanfhyp = (`meanthetabl',`meanthetabl'+`betanfhyp_theta',`meanfatigbl', `meanfatigbl'+`betanfhyp_fatig')
	    
        * create a vector that contains the standard deviations
	    matrix sdthetanfhyp = (`sdthetabl',`sdthetabl',`sdfatigbl',`sdfatigbl')
	    * draw a sample of specified cases from a normal distribution with specified correlation structure
	    * and specified means and standard deviations
	    drawnorm thetabl thetapost fatigbl fatigpost, n(`nperg') corr(cthetanfhyp) cstorage(lower) means(mthetanfhyp) sds(sdthetanfhyp) double
        gen tx=1
        save datanfhyp, replace
 
        drop _all
        use datawl
        append using datanfhyp
        save datatemp.dta,replace
 */       
        
        gen id=_n
        
        
        **** end of generating data

        *** begin analyzing data

        count
        local totaln=r(N)

        sum  deltaY if tx==0
        local meandeltaY0=r(mean)
        local sddeltaY0=r(sd)
        sum  deltaY if tx==1
        local meandeltaY1=r(mean)
        local sddeltaY1=r(sd)

        sum  M1 if tx==0
        local meanM10=r(mean)
        local sdM10=r(sd)
        sum  M1 if tx==1
        local meanM11=r(mean)
        local sdM11=r(sd)

        sum  M2 if tx==0
        local meanM20=r(mean)
        local sdM20=r(sd)
        sum  M2 if tx==1
        local meanM21=r(mean)
        local sdM21=r(sd)
        
        sum  M3 if tx==0
        local meanM30=r(mean)
        local sdM30=r(sd)
        sum  M3 if tx==1
        local meanM31=r(mean)
        local sdM31=r(sd)
        
        sum  M4 if tx==0
        local meanM40=r(mean)
        local sdM40=r(sd)
        sum  M4 if tx==1
        local meanM41=r(mean)
        local sdM41=r(sd)
        
        sum  M5 if tx==0
        local meanM50=r(mean)
        local sdM50=r(sd)
        sum  M5 if tx==1
        local meanM51=r(mean)
        local sdM51=r(sd)
        
        sum tx 
        local meantx=r(mean)
        local sdtx=r(sd)
      
        *** models for mediation of tx on change in pain intensity (deltaY) by M1 - M5

        sureg (M1 tx) (M2 tx) (M3 tx) (M4 tx) (M5 tx) (deltaY M1 M2 M3 M4 M5 tx)
        matrix list r(table)
        local a1=r(table)[1,1]
        local a2=r(table)[1,3]
        local a3=r(table)[1,5]
        local a4=r(table)[1,7]
        local a5=r(table)[1,9]
        local b1=r(table)[1,11]
        local b2=r(table)[1,12]
        local b3=r(table)[1,13]
        local b4=r(table)[1,14]
        local b5=r(table)[1,15]
        local cprime=r(table)[1,16]

        nlcom [M1]_b[tx]*[deltaY]_b[M1]
        matrix list r(b)
        matrix list r(V)
        display sqrt(r(V)[1,1])
        local a1b1=r(b)[1,1]
        local a1b1var=r(V)[1,1]
        local a1b1sd=sqrt(`a1b1var')
        local a1b1upp=`a1b1'+1.959964*`a1b1sd'
        local a1b1no0=(`a1b1upp'<0)
        display `a1b1no0'
        
        nlcom [M2]_b[tx]*[deltaY]_b[M2]
        matrix list r(b)
        matrix list r(V)
        display sqrt(r(V)[1,1])
        local a2b2=r(b)[1,1]
        local a2b2var=r(V)[1,1]
        local a2b2sd=sqrt(`a2b2var')
        local a2b2upp=`a2b2'+1.959964*`a2b2sd'
        local a2b2no0=(`a2b2upp'<0)

        nlcom [M3]_b[tx]*[deltaY]_b[M3]
        matrix list r(b)
        matrix list r(V)
        display sqrt(r(V)[1,1])
        local a3b3=r(b)[1,1]
        local a3b3var=r(V)[1,1]
        local a3b3sd=sqrt(`a3b3var')
        local a3b3upp=`a3b3'+1.959964*`a3b3sd'
        local a3b3no0=(`a3b3upp'<0)

        nlcom [M4]_b[tx]*[deltaY]_b[M4]
        matrix list r(b)
        matrix list r(V)
        display sqrt(r(V)[1,1])
        local a4b4=r(b)[1,1]
        local a4b4var=r(V)[1,1]
        local a4b4sd=sqrt(`a4b4var')
        local a4b4upp=`a4b4'+1.959964*`a4b4sd'
        local a4b4no0=(`a4b4upp'<0)

        nlcom [M5]_b[tx]*[deltaY]_b[M5]
        matrix list r(b)
        matrix list r(V)
        display sqrt(r(V)[1,1])
        local a5b5=r(b)[1,1]
        local a5b5var=r(V)[1,1]
        local a5b5sd=sqrt(`a5b5var')
        local a5b5upp=`a5b5'+1.959964*`a5b5sd'
        local a5b5no0=(`a5b5upp'<0)

        nlcom [M1]_b[tx]*[deltaY]_b[M1]+[M2]_b[tx]*[deltaY]_b[M2]+[M3]_b[tx]*[deltaY]_b[M3]+[M4]_b[tx]*[deltaY]_b[M4]+[M5]_b[tx]*[deltaY]_b[M5]
        matrix list r(b)
        matrix list r(V)
        display sqrt(r(V)[1,1])
        local ab=r(b)[1,1]
        local abvar=r(V)[1,1]
        local absd=sqrt(`abvar')
        local abupp=`ab'+1.959964*`absd'
        local abno0=(`abupp'<0)
        display `abno0'

        regress deltaY tx
        matrix list r(table)
        local toteffmodel3=r(table)[1,1]
        display `toteffmodel3'

        local percind=`ab'/`toteffmodel3'
        display `percind'  

*        bootstrap r(indM1) r(indM2) r(indM3) r(indM4) r(indM5) r(indtotal), saving(bootsav, replace) bca reps(1000): bootmpm
*        return list
*        matrix list r(table)

*        estat boot, percentile bc bca
*        return list


*       *** for testing
*        if `nummerin'==44 {
*         capture save data`nummerin', replace
*        }
     
        return scalar totaln=`totaln'

        return scalar meandeltaY0=`meandeltaY0'
        return scalar sddeltaY0=`sddeltaY0'
        return scalar meandeltaY1=`meandeltaY1'
        return scalar sddeltaY1=`sddeltaY1'

        return scalar meanM10=`meanM10'
        return scalar sdM10=`sdM10'
        return scalar meanM11=`meanM11'
        return scalar sdM11=`sdM11'

        return scalar meanM20=`meanM20'
        return scalar sdM20=`sdM20'
        return scalar meanM21=`meanM21'
        return scalar sdM21=`sdM21'

        return scalar meanM30=`meanM30'
        return scalar sdM30=`sdM30'
        return scalar meanM31=`meanM31'
        return scalar sdM31=`sdM31'

        return scalar meanM40=`meanM40'
        return scalar sdM40=`sdM40'
        return scalar meanM41=`meanM41'
        return scalar sdM41=`sdM41'

        return scalar meanM50=`meanM50'
        return scalar sdM50=`sdM50'
        return scalar meanM51=`meanM51'
        return scalar sdM51=`sdM51'
        
        return scalar meantx=`meantx'
        return scalar sdtx=`sdtx'
        
        return scalar a1=`a1'
        return scalar a2=`a2'
        return scalar a3=`a3'
        return scalar a4=`a4'
        return scalar a5=`a5'
        return scalar b1=`b1'
        return scalar b2=`b2'
        return scalar b3=`b3'
        return scalar b4=`b4'
        return scalar b5=`b5'
        return scalar cprime=`cprime'

        return scalar a1b1=`a1b1'
        return scalar a1b1sd=`a1b1sd'
        return scalar a1b1upp=`a1b1upp'
        return scalar a1b1no0=`a1b1no0'
                
        return scalar a2b2=`a2b2'
        return scalar a2b2sd=`a2b2sd'
        return scalar a2b2upp=`a2b2upp'
        return scalar a2b2no0=`a2b2no0'

        return scalar a3b3=`a3b3'
        return scalar a3b3sd=`a3b3sd'
        return scalar a3b3upp=`a3b3upp'
        return scalar a3b3no0=`a3b3no0'

        return scalar a4b4=`a4b4'
        return scalar a4b4sd=`a4b4sd'
        return scalar a4b4upp=`a4b4upp'
        return scalar a4b4no0=`a4b4no0'

        return scalar a5b5=`a5b5'
        return scalar a5b5sd=`a5b5sd'
        return scalar a5b5upp=`a5b5upp'
        return scalar a5b5no0=`a5b5no0'

        return scalar ab=`ab'
        return scalar absd=`absd'
        return scalar abupp=`abupp'
        return scalar abno0=`abno0'

        return scalar toteffmodel3=`toteffmodel3'
        return scalar percind=`percind'
        
        save tempdata, replace

    end        

* see mediation5-shell for how this program is called


     
        
