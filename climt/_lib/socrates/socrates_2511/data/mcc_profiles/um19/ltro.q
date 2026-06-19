netcdf ltro.q{                                                                  

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
                q =   .325312E-05,  .324904E-05,  .325225E-05,  .324747E-05,
                      .326490E-05,  .389924E-05,  .135873E-04,  .565319E-04,
                      .167071E-03,  .382036E-03,  .774811E-03,  .152610E-02,
                      .287933E-02,  .520183E-02,  .945657E-02,  .115333E-01,
                      .135005E-01,  .152113E-01,  .161483E-01;

}                                                                               
