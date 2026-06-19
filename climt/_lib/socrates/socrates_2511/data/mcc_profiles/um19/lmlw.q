netcdf lmlw.q{                                                                  

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
             plev =   .534450E+03,  .152700E+04,  .305400E+04,  .585350E+04,
                      .101800E+05,  .152700E+05,  .203600E+05,  .254500E+05,
                      .305400E+05,  .361390E+05,  .430105E+05,  .514090E+05,
                      .610800E+05,  .712600E+05,  .806765E+05,  .885659E+05,
                      .947248E+05,  .992550E+05,  .101495E+06;
                q =   .400013E-05,  .399912E-05,  .400206E-05,  .400361E-05,
                      .408790E-05,  .606046E-05,  .106520E-04,  .182665E-04,
                      .380812E-04,  .839007E-04,  .203190E-03,  .444080E-03,
                      .844138E-03,  .133649E-02,  .179157E-02,  .209851E-02,
                      .234999E-02,  .255933E-02,  .267399E-02;

}                                                                               
