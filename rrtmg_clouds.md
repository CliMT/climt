# RRMTG Clouds with McICA

## Cloud properties

There are three options for the RRTMG `inflag`, as given in the `climt` dictionary:

```python
rrtmg_cloud_props_dict = {
    'direct_input': 0,
    'single_cloud_type': 1,
    'liquid_and_ice_clouds': 2
}
```

With McICA, we cannot use `single_cloud_type`, but can choose between `direct_input` and `liquid_and_ice_clouds`.
If we choose `direct_input`, we input the `longwave_optical_thickness_due_to_cloud`, `shortwave_optical_thickness_due_to_cloud`, as well as the shortwave parameters `single_scattering_albedo_due_to_cloud`, `cloud_asymmetry_parameter` and `cloud_forward_scattering_fraction`.
The `cloud_forward_scattering_fraction` is used to scale the other shortwave parameters (`shortwave_optical_thickness_due_to_cloud`, `single_scattering_albedo_due_to_cloud` and `cloud_asymmetry_parameter`), but it is not directly used in the radiative transfer calculations.
The other cloud properties, namely `cloud_ice_particle_size` and `cloud_water_droplet_radius`, `mass_content_of_cloud_ice_in_atmosphere_layer`, and `mass_content_of_cloud_liquid_water_in_atmosphere_layer` are completely unused.
The RRTMG `iceflag` and `liqflag` are irrelevant.

On the other hand, if we choose `liquid_and_ice_clouds`,
any input values for `longwave_optical_thickness_due_to_cloud`, `shortwave_optical_thickness_due_to_cloud`, `single_scattering_albedo_due_to_cloud` and `cloud_asymmetry_parameter` are irrelevant.
Instead, these parameters are calculated from the cloud ice and water droplet particle sizes (`cloud_ice_particle_size` and `cloud_water_droplet_radius`), as well as the cloud ice and water paths (`mass_content_of_cloud_ice_in_atmosphere_layer`, `mass_content_of_cloud_liquid_water_in_atmosphere_layer`).
The methods used for the calculations depend on the cloud ice and water properties; `iceflag` and `liqflag` in RRTMG.


## Calculation of cloud properties

### Longwave

For the longwave, the only cloud property of interest for calculating radiative fluxes in RRTMG, is the optical depth. This is calculated as follows.

```fortran
do lay = 1, nlayers
  do ig = 1, ngptlw
    taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) + &
                     clwpmc(ig,lay) * abscoliq(ig)
  end do
end do
```

Values of cloud optical depth `taucmc` are calculated for each model layer (pressure), `lay`, and each g-interval, `ig`.
The cloud ice and liquid absorption coefficients (`abscoice` and `abscoliq`) are multiplied by the cloud ice and liquid water paths (`ciwpmc` and `clwpmc`) respectively, to give the ice cloud optical depth and the liquid water cloud optical depth.
The cloud ice and liquid water paths are input by the user, in `climt` as `mass_content_of_cloud_ice_in_atmosphere_layer` and `mass_content_of_cloud_liquid_water_in_atmosphere_layer` respectively.
The cloud ice and liquid absorption coefficients are calculated based on the ice and liquid water particle sizes (specified by the user), and this calculation depends on the choice of `iceflag` and `liqflag`.

### Shortwave

For the shortwave, there are three cloud properties, which affect the radiative flux calculation in RRTMG, namely the optical depth, the single scattering albedo and the asymmetry parameter.

The shortwave optical depth is calculated as:
```fortran
taucmc(ig,lay) = tauice + tauliq
```
with
```fortran
tauice = (1 - forwice(ig) + ssacoice(ig)) * ciwpmc(ig,lay) * extcoice(ig)
tauliq = (1 - forwliq(ig) + ssacoliq(ig)) * clwpmc(ig,lay) * extcoliq(ig)
```

The single scattering albedo is calculated as:
```fortran
ssacmc(ig,lay) = (scatice + scatliq) / taucmc(ig,lay)
```
with
```fortran
scatice = ssacoice(ig) * (1._rb - forwice(ig)) / (1._rb - forwice(ig) * ssacoice(ig)) * tauice
scatliq = ssacoliq(ig) * (1._rb - forwliq(ig)) / (1._rb - forwliq(ig) * ssacoliq(ig)) * tauliq
```

The asymmetry parameter is given by:
```fortran
asmcmc(ig,lay) = 1.0_rb / (scatliq + scatice) * ( &
                   scatliq * (gliq(ig) - forwliq(ig)) / (1.0_rb - forwliq(ig)) + &
                   scatice * (gice(ig) - forwice(ig)) / (1.0_rb - forwice(ig)) )
```

Values of optical depth, single scattering albedo and asymmetry parameter are calculated for each model layer (pressure), `lay`, and each g-interval, `ig`. The cloud ice and liquid water paths (`ciwpmc` and `clwpmc`) are input by the user.
The other parameters (`extcoice`, `extcoliq`, `ssacoice`, `ssacoliq`, `gice`, `gliq`, `forwice`, `forwliq`) are calculated based on the ice and liquid water particle sizes and this calculation depends on the choice of `iceflag` and `liqflag`.

## Cloud ice properties

There are four options for the RRTMG `iceflag`. These are given in the `climt` dictionary:

```python
rrtmg_cloud_ice_props_dict = {
    'ebert_curry_one': 0,
    'ebert_curry_two': 1,
    'key_streamer_manual': 2,
    'fu': 3
}
```

### ebert_curry_one

For the longwave, `ebert_curry_one` gives an absorption coefficient of
```
abscoice = 0.005 + 1.0 / radice
```
Here, `radice` is the ice particle size and the absorption coefficient is the same for all wavebands.

`ebert_curry_one` should not be used for the shortwave component with McICA.

### ebert_curry_two

`ebert_curry_two` is the default choice for cloud ice optical properties in `climt`.
In this case, the longwave absorption coefficient is calculated as follows.
```fortran
do ig = 1, ngptlw
  ib = icb(ngb(ig))
  abscoice(ig) = absice1(1,ib) + absice1(2,ib)/radice
```
The absorption coefficient `abscoice` is a function of g-interval `ig` and is made up of two contributions.
The first of these `absice1(1, ib)` comes from a look up table and is given in [m<sup>2</sup> g].
`ib` provides an index for the look up table, based on the waveband of the g-interval.
`absice1(2,ib)` also comes from a look up table and is given in [microns m<sup>2</sup> g]. It is divided by `radice`, the cloud ice particle size, providing an ice particle size dependence of the absorption coefficient.
Although the syntax does not emphasise it, the absorption coefficient may also depend on model layer (pressure), as `radice` can have model layer dependence.
`radice` comes from the input property labeled `cloud_ice_particle_size` in `climt`.
The ice particle size dependent term is more important than the independent term (`absice1(2,ib)/radice` > `absice1(1,ib)`) at all wavebands for ice particle sizes less than 88 microns.
Using `ebert_curry_two`, the ice particle size must be in the range [13, 130] microns, and even for larger particle sizes (> 88), the ice particle size dependent term is more important than the independent term for four of the five wavebands.

For the shortwave, the parameters (required for the optical depth, single scattering albedo and asymmetry) are calculated as follows.
```fortran
extcoice(ig) = abari(icx) + bbari(icx)/radice
ssacoice(ig) = 1._rb - cbari(icx) - dbari(icx) * radice
gice(ig) = ebari(icx) + fbari(icx) * radice
forwice(ig) = gice(ig)*gice(ig)
```
`abari`, `bbari`, `cbari`, `dbari`, `ebari` and `fbari` are all look up tables, containing five values, which correspond to five different wavebands.
The choice of waveband is indicated by `icx`.
The particle size dependence comes from `radice`, so each parameter consists of both a size independent and a size dependent contribution.


### key_streamer_manual

In this case, both the longwave absorption coefficient and three of the shortwave parameters (`excoice`, `ssacoice`, `gice`) are interpolated from look up tables.
Comments in the RRTMG code state that these look up tables are for a spherical ice particle parameterisation.
The look up tables contain 42 values for each of the 16 longwave and 14 shortwave wavebands.
The 42 values correspond to different ice particle radii, evenly spaced in the range [5, 131] microns.
Ice particles must be within this range, otherwise an error is thrown.

The shortwave parameter `forwice` is calculated as the square of `gice`.

### fu

The longwave absorption coefficient and shortwave parameters `extcoice`, `ssacoice` and `gice` are interpolated from look up tables.
The look up tables differ to those in `key_streamer_manual`, and
comments in the RRTMG code state that the look up tables for `fu` are for a hexagonal ice particle parameterisation.
The look up tables for `fu` are slightly larger than those for `key_streamer_manual`, and the range of allowed values for the ice particle size is corresponding larger ([5, 140] microns).

The shortwave parameter `forwice` is calculated from `fdelta` (again taken from look up tables) and `ssacoice` as follows.

```
forwice = fdelta + 0.5 / ssacoice
```

## Cloud liquid properties

There are two options for the RRTMG `liqflag`. These are given in the `climt` dictionary:

```python
rrtmg_cloud_liquid_props_dict = {
    'radius_independent_absorption': 0,
    'radius_dependent_absorption': 1
}
```

For `radius_independent_absorption`, the longwave absorption coefficient is 0.0903614 for all wavebands.
This option should not be used for the shortwave.

`radius_dependent_absorption` is the default choice for cloud liquid water properties in `climt`.
In this case, the longwave absorption coefficient and the shortwave parameters `extcoliq`, `ssacoliq`, and `gliq` are interpolated from look up tables.
The look up tables have values for particle sizes in the range [2.5, 59.5] microns in 1 micron intervals (58 values) for each of the 16 longwave and 14 shortwave wavebands.
The shortwave parameter `forwliq` is calculated as the square of `gliq`.

## Cloud overlap method

This is the RRTMG `icld`.

```python
rrtmg_cloud_overlap_method_dict = {
    'clear_only': 0,
    'random': 1,
    'maximum_random': 2,
    'maximum': 3
}
```

If we choose `clear_only`, there are no clouds, regardless of the other input.
There is no difference between the other three options (`random`, `maximum_random` and `maximum`) for the McICA version of RRTMG.
