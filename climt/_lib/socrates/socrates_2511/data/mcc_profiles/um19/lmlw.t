netcdf lmlw.t{                                                                  

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
             plev =   .534450E+03,  .152700E+04,  .305400E+04,  .585350E+04,
                      .101800E+05,  .152700E+05,  .203600E+05,  .254500E+05,
                      .305400E+05,  .361390E+05,  .430105E+05,  .514090E+05,
                      .610800E+05,  .712600E+05,  .806765E+05,  .885659E+05,
                      .947248E+05,  .992550E+05,  .101495E+06;
                t =   .219529E+03,  .215777E+03,  .215200E+03,  .216043E+03,
                      .217400E+03,  .218480E+03,  .219293E+03,  .223048E+03,
                      .229714E+03,  .236821E+03,  .244371E+03,  .252100E+03,
                      .259223E+03,  .264139E+03,  .267087E+03,  .269283E+03,
                      .270852E+03,  .271807E+03,  .272158E+03;

}                                                                               
