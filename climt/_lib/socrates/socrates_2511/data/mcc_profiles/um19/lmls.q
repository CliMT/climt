netcdf lmls.q{                                                                  

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  19;


variables:
    float lat(lat);                                                           
             lat:units = "degree";                                            
             lat:title = "LATITUDE";                                         
    float lon(lon);                                                           
             lon:units = "degree";                                            
             lon:title = "LONGITUDE";                                        
    float plev(plev);                                                         
             plev:units = "Pa";                                               
             plev:title = "PRESSURE";                                        

    float q(plev,lon,lat);                                                     
             q:units = "None";                                                
             q:title = "MMR OF WATER VAPOUR";                                                

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .531825E+03,  .151950E+04,  .303900E+04,  .582475E+04,
                      .101300E+05,  .151950E+05,  .202600E+05,  .253250E+05,
                      .303900E+05,  .359615E+05,  .427993E+05,  .511565E+05,
                      .607800E+05,  .709100E+05,  .802802E+05,  .881309E+05,
                      .942596E+05,  .987675E+05,  .100996E+06;
                q =   .399661E-05,  .400214E-05,  .400981E-05,  .399821E-05,
                      .400173E-05,  .428756E-05,  .163013E-04,  .796482E-04,
                      .201571E-03,  .357356E-03,  .645038E-03,  .108591E-02,
                      .205899E-02,  .385590E-02,  .601352E-02,  .799227E-02,
                      .963579E-02,  .109396E-01,  .116537E-01;

}                                                                               
