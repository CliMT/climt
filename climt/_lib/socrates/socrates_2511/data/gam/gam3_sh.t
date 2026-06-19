netcdf gam3_sh.t{                                                               

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  18;


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
             plev =   .100000E+03,  .300000E+03,  .100000E+04,  .300000E+04,
                      .500000E+04,  .700000E+04,  .100000E+05,  .150000E+05,
                      .200000E+05,  .203471E+05,  .205516E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .983021E+05;
                t =   .265430E+03,  .252160E+03,  .225330E+03,  .218680E+03,
                      .215700E+03,  .214130E+03,  .213980E+03,  .215730E+03,
                      .217600E+03,  .217790E+03,  .217900E+03,  .220060E+03,
                      .225480E+03,  .238150E+03,  .249720E+03,  .265260E+03,
                      .273590E+03,  .276570E+03;

}                                                                               
