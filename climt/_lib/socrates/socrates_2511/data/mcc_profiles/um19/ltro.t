netcdf ltro.t{                                                                  

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
                t =   .233001E+03,  .223294E+03,  .211958E+03,  .200404E+03,
                      .202936E+03,  .215444E+03,  .226316E+03,  .235433E+03,
                      .243777E+03,  .252097E+03,  .261168E+03,  .270317E+03,
                      .279202E+03,  .285754E+03,  .290343E+03,  .294460E+03,
                      .297435E+03,  .299252E+03,  .299920E+03;

}                                                                               
