netcdf gam3_nh.t{                                                               

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
                      .200000E+05,  .202088E+05,  .204119E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .944467E+05;
                t =   .264420E+03,  .251260E+03,  .226220E+03,  .220140E+03,
                      .216970E+03,  .215700E+03,  .215370E+03,  .217600E+03,
                      .219500E+03,  .219630E+03,  .219740E+03,  .222110E+03,
                      .227680E+03,  .240640E+03,  .252000E+03,  .267250E+03,
                      .274670E+03,  .277570E+03;

}                                                                               
