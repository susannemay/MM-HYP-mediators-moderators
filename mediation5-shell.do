
**************************************************************************
*
* Program to try out power and sample size for project with Mark Jensen and Melissa Day
*
* this is the "shell" program for mediation5.do 
*
*
* Author: Susanne May
* Date: 09/19/2020 
*
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
**************************************************************************
    display "Simulation started  $S_DATE   $S_TIME"
    
    clear
    use one
    replace nummer=1 in 1
    save,replace
    clear

    simulate     totaln=r(totaln)  ///
        meandeltaY0=r(meandeltaY0) ///
        sddeltaY0=r(sddeltaY0) ///
        meandeltaY1=r(meandeltaY1) ///
        sddeltaY1=r(sddeltaY1) ///
        meanM10=r(meanM10)  ///
        sdM10=r(sdM10)  ///
        meanM11=r(meanM11)  ///
        sdM11=r(sdM11)  ///
        meanM20=r(meanM20)  ///
        sdM20=r(sdM20)  ///
        meanM21=r(meanM21)  ///
        sdM21=r(sdM21)  ///
        meanM30=r(meanM30)  ///
        sdM30=r(sdM30)  ///
        meanM31=r(meanM31)  ///
        sdM31=r(sdM31)  ///
        meanM40=r(meanM40)  ///
        sdM40=r(sdM40)  ///
        meanM41=r(meanM41)  ///
        sdM41=r(sdM41)  ///
        meanM50=r(meanM50)  ///
        sdM50=r(sdM50)  ///
        meanM51=r(meanM51)  ///
        sdM51=r(sdM51)  ///
        meantx=r(meantx)  ///
        sdtx=r(sdtx)  ///
        a1=r(a1)  ///
        a2=r(a2)  ///
        a3=r(a3)  ///
        a4=r(a4)  ///
        a5=r(a5)  ///
        b1=r(b1)  ///
        b2=r(b2)  ///
        b3=r(b3)  ///
        b4=r(b4)  ///
        b5=r(b5)  ///
        cprime=r(cprime)  ///
        a1b1=r(a1b1)  ///
        a1b1sd=r(a1b1sd)  ///
        a1b1upp=r(a1b1upp)  ///
        a1b1no0=r(a1b1no0)  ///
        a2b2=r(a2b2)  ///
        a2b2sd=r(a2b2sd)  ///
        a2b2upp=r(a2b2upp)  ///
        a2b2no0=r(a2b2no0)  ///
        a3b3=r(a3b3)  ///
        a3b3sd=r(a3b3sd)  ///
        a3b3upp=r(a3b3upp)  ///
        a3b3no0=r(a3b3no0)  ///
        a4b4=r(a4b4)  ///
        a4b4sd=r(a4b4sd)  ///
        a4b4upp=r(a4b4upp)  ///
        a4b4no0=r(a4b4no0)  ///
        a5b5=r(a5b5)  ///
        a5b5sd=r(a5b5sd)  ///
        a5b5upp=r(a5b5upp)  ///
        a5b5no0=r(a5b5no0)  ///
        ab=r(ab)  ///
        absd=r(absd)  ///
        abupp=r(abupp)  ///
        abno0=r(abno0)  ///
        toteffmodel3=r(toteffmodel3)  ///
        percind=r(percind)  ///
	    , reps(`6'): simula, nperg(`1') toteff(`2') indpoftot(`3') toteffsd(`4') seed(`5')

***
    display "Program was run as: do mediation5-shell `*'"
    display "with arguments: nperg toteff indpoftot toteffsd seed"
    sum 
    gen powerab=abno0
    
    list powerab abno0 percind  ab cprime toteffmodel3 if percind>1
        
    sum powerab
    sum powerab if percind<1
    

    capture save results.dta, replace

    display "Simulation ended  $S_DATE   $S_TIME"
    
    



