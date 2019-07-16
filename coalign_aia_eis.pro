;#NAME: coalign_aia_eis
;
;#PURPOSE: to coalign data from SDO/AIA and Hinode/EIS using 'eis_aia_offset.pro'
;
;#CALLING SEQUENCE: coalign_aia_eis, aia_cube, eis_dir, wave=wave
;
;#KEYWORDS: aia_cube: an AIA datacube containing image sequences (full Sun FOV)
;           eis_dir: the directory containing EIS datacube for a specific spectral line
;           wave: the passband of AIA, in string type, e.g. '171' for AIA 171 passband.
;
;#OUTPUT: save files containing the shifted AIA data, observing time and the EIS field-of-view parameters.
;
;#CAUTION: Here the code is specified for the datasets I used, in order to make it flexible, you should
;           change the x and y pixel numbers and EIS observing time, i.e. solar_x[79], solar_y[303] and time[40]
;            to your own data parameters
;
;#AUTHOR: Zihao Yang, Peking Univerisity, July 11th, 2019, at UCL-MSSL.

pro coalign_aia_eis, aia_cube, eis_dir, wave=wave

restore, aia_cube
aia_index=index
aia_int=int 

f_eis=file_search(eis_dir+'./sgf*.sav')
fc=n_elements(f_eis)

for jj=0, fc-1 do begin
  restore, f_eis[jj]
  xc_eis=(solar_x[0]+solar_x[79])/2
  yc_eis=(solar_y[0]+solar_y[303])/2
  t_eis=time[40]
  dx_eis=(solar_x[79]-solar_x[0])/79
  dy_eis=(solar_y[303]-solar_y[0])/303


  aia_i=aia_index[jj]
  aia_d=aia_int[*,*,jj]
  aia_time=aia_i.date_obs

  nx=aia_i.naxis1
  ny=aia_i.naxis2
  dx=aia_i.cdelt1
  dy=aia_i.cdelt2

  solar_x_aia=findgen(nx)*dx-nx/2*dx
  solar_y_aia=findgen(ny)*dy-ny/2*dy

  offsets=eis_aia_offsets(t_eis)
  aia_shift=shift(aia_d, -(offsets[0]*dx_eis)+10, -(offsets[1]*dy_eis)-3)

  x_pix=fltarr(2)
  x_pix[0]=(solar_x[0]-solar_x_aia[0])/dx
  x_pix[1]=(solar_x[79]-solar_x_aia[0])/dx

  y_pix=fltarr(2)
  y_pix[0]=(solar_y[0]-solar_y_aia[0])/dy
  y_pix[1]=(solar_y[303]-solar_y_aia[0])/dy

  save,filename='./'+wave+'/AIA_'+wave+'_coaligned'+strtrim(jj+10,2)+'.sav',aia_shift,x_pix,y_pix,aia_time
endfor
end