netcdf lmls.t{                                                                  

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

    float t(plev,lon,lat);                                                     
             t:units = "K";                                                   
             t:title = "TEMPERATURE";                                                        

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .531825E+03,  .151950E+04,  .303900E+04,  .582475E+04,
                      .101300E+05,  .151950E+05,  .202600E+05,  .253250E+05,
                      .303900E+05,  .359615E+05,  .427993E+05,  .511565E+05,
                      .607800E+05,  .709100E+05,  .802802E+05,  .881309E+05,
                      .942596E+05,  .987675E+05,  .100996E+06;
                t =   .236562E+03,  .227242E+03,  .220479E+03,  .216969E+03,
                      .216000E+03,  .218366E+03,  .225662E+03,  .234608E+03,
                      .242528E+03,  .250813E+03,  .259239E+03,  .267404E+03,
                      .275262E+03,  .282065E+03,  .287064E+03,  .290278E+03,
                      .292325E+03,  .293512E+03,  .293948E+03;

}                                                                               
