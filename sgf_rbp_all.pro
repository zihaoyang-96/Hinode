; NAME : sgf_rbp_all
; 
; PURPOSE : Perform single Gaussian fit and derive RBp profile (a modified version of the RB asymmetry analysis originally developped by 
;                    De Pontieu et al. 2009, ApJ, 701, L1; see the definition in Section 2 of Tian et al. 2011, ApJ, 738, 18) for line profiles
;
; CALLING SEQUENCE : sgf_rbp_all
; 
;      
; KEYWORDS : wave0 - rest wavelength, use the average centroid if not set 
;         
;            posi - position of the line profile plot
;                    
;            xtitle - lable of x-cooridinate of the line profile plot and RB profile plot
;                    
;            xtitleRB -  lable of y-cooridinate of the RB profile plot
;                    
;            ytitle -  lable of y-cooridinate of the line profile plot
;                    
;            xav - running average of the profiles over xav*2+1 pixels in x direction, no avearge if xav=0 or not set
;                    
;            yav - running average of the profiles over yav*2+1 pixels in y direction, no avearge if yav=0 or not set
; 
; OUTPUTS : fit_para - 3-D array[3,nx,ny] of single-Gaussian fit parameters: a[0], a[1], a[2]

;          para_sgf -  3-D array[4,nx,ny] of line parameters: peak intensity, Doppler shift, line width, background intensity
;          
;          chisq_sgf - 2-D array[nx,ny] of the chisq value for SGF
;                    
;          para_rb - 3-D array[3,nx,ny] of RBp parameters: relative intensity, velocity, and width
;                    
;          rb_array - RBp profiles [nx,ny,15]
;                    
;          r_array - relative red wing emission as a fuction of velocity [nx,ny,15]
;                    
;          b_array - relative blue wing emission as a fuction of velocity [nx,ny,15]
;
; HISTORY : First edition by Hui Tian at CfA, April 2, 2013
;          Modified by Zihao Yang (PKU) at UCL-MSSL, July 16, 2019 (add the fit_para and modify some parts to make it robust)
;
; NOTES: the data cubes read by sgf_rbp_all.pro should contain a 3D array of intensity for each pixel [nw,nx,ny] (nw is the number of wavelength locations)

;
;; **********************************************************************

pro sgf_rbp_all,wave0=wave0,posi=posi,xtitle=xtitle,xtitleRB=xtitleRB,ytitle=ytitle,xav=xav,yav=yav

!p.font=-1
f=file_search('./*datacube*.sav', count = fc)

if not(keyword_set(xtitle)) then xtitle = 'Wavelength (!N!6!sA!r!u!9 %!6 !n!N!3)'
if not(keyword_set(xtitleRB)) then xtitleRB = 'Velocity (km/s)'
if not(keyword_set(ytitle)) then ytitle = 'Counts';'Radiance (erg !Ns!E-1!Ncm!E-2!Nsr!E-1 !N!6!sA!r!u!9 %!6!n!U-1!N!3)'
if not(keyword_set(posi)) then posi=[0.07,0.10,0.60,0.95]

for num=0,fc-1 do begin
file_lps=f[num]
restore,file_lps
dlambda=mean(deriv(wvl)) ;spectral resolution
nw=(size(int))[1] 
nx=(size(int))[2] 
ny=(size(int))[3]

bad=where(int le 0)
inten=int
inten[bad]='nan'
;if not(keyword_set(wave0)) then begin
; centro=fltarr(11,21)
; for xx=0,10 do begin
; for yy=0,20 do begin
;   res=mpfitpeak(wvl,inten[*,xx+5,yy+270], a, nterms = 4, /double, /positive) 
;   ; res = mpfitpeak(wvl,average(average(inten[*,30:40,20:30],2),2,missing=NaN), a, nterms = 4, /double, /positive) 
;   ;res = mpfitpeak(wvl,average(average(inten,2),2), a, nterms = 4, /double, /positive) 
;   ;res = mpfitpeak(wvl,average(average(inten[*,5:15,270:300],2),2,missing=NaN), a, nterms = 4, /double, /positive) 
res = mpfitpeak(wvl,average(average(inten,2),2,missing=NaN), a, nterms = 4, /double, /positive) 
;   centro[xx,yy]= a[1] ;set the rest wavelength by assuming 0 shift of the average profile
;   endfor
;   endfor
wave0=a[1]
;wave0=average(centro)
;endif
print,wave0
; stop

loadct,0  & tvlct,255L,0L,0L,2L  & tvlct,0L,255L,0L,3L   & tvlct,0L,0L,255L,4L
window, 1, xs =1200, ys =600, xpos=20, ypos=100

para_sgf=fltarr(4,nx,ny)
chisq_sgf=fltarr(nx,ny)
fit_para=fltarr(3,nx,ny)

steps = findgen(15)*10. ; R-B Steps in km/s
rb_array=fltarr(nx,ny,15)
r_array=fltarr(nx,ny,15)
b_array=fltarr(nx,ny,15)
para_rb=fltarr(3,nx,ny)

;average of the spectra over (xav*2+1)*(yav*2+1) pixels
if not(keyword_set(xav)) then xav=0
if not(keyword_set(yav)) then yav=0

for m=0, ny-1 do begin 
for n=0, nx-1 do begin
if (n eq -1) then begin
;if (n eq 46) or (n eq 57) or (n eq 66) or (n eq 67) or (n eq 68) or (n eq 78) or (n eq 79) then begin
para_sgf[0,n,m]=0
para_sgf[1,n,m]=0
para_sgf[2,n,m]=0
para_sgf[3,n,m]=0
chisq_sgf[n,m]=0
para_rb[0,n,m]=0
para_rb[1,n,m]=0
para_rb[2,n,m]=0
rb_array[n,m,*]=0
fit_para[*,n,m]=0
endif else begin
;if no bad stripes then comment this line
y=dblarr(nw)  &  err_ave=dblarr(nw) 
x1=(n-xav) > 0 
x2=(n+xav) < (nx-1) 
y1=(m-yav) > 0 
y2=(m+yav) < (ny-1) 
for iw=0,nw-1 do begin
  if (x1 eq x2) and (y1 eq y2) then begin
    y[iw]=reform(int[iw,x1,y1])  
    err_ave[iw]=reform(err[iw,x1,y1])
    goto,nextprofile
  endif
  if (x1 eq x2) and (y1 ne y2) then begin
    int_tmp=int[iw,x1,y1:y2]  
    err_tmp=err[iw,x1,y1:y2]
  endif
  if (x1 ne x2) and (y1 eq y2) then begin
    int_tmp=int[iw,x1:x2,y1]  
    err_tmp=err[iw,x1:x2,y1]
  endif
  if (x1 ne x2) and (y1 ne y2) then begin
    int_tmp=int[iw,x1:x2,*]  &  int_tmp=int_tmp[*,*,y1:y2]
    err_tmp=err[iw,x1:x2,*]  &  err_tmp=err_tmp[*,*,y1:y2]
  endif
  ;not to avearge bad data
  sub_good=where(int_tmp ne -32768,num_good)
  ;sub_good=where(err_tmp ne -100,num_good)
  if num_good gt 0 then begin
    y[iw]=average(int_tmp[sub_good])
    err_ave[iw]=sqrt(total(err_tmp[sub_good]^2))/num_good
  endif
  nextprofile:
endfor
err_ave=err_ave>0.001<(abs(y)*0.999) ;make sure that the error value is larger than 0 and smaller than the data value

;single Gaussian fit
x=wvl
plot_io,x,y,title='SGF to '+lineid+' line profile',position=posi,xtitle=xtitle,ytitle=ytitle,psym=4,$
      xrange=[min(x)-dlambda/2,max(x)+dlambda/2],xstyle=1,yrange=[min(y)*0.95>0.1,max(y)*1.5],ystyle=1,charsize=1.5
err_plot,x,y-err_ave,y+err_ave,width=0.005,thick=1
;bpp = where(y gt 0., bppc)  &  ee = err_ave
bpp = where(y gt (-50), bppc)  &  ee = err_ave ;sometime the dark current subtraction is not perfect
if bppc le 4 then begin ;do not do fit
  a=fltarr(4)-999. & chisq=-999. & i_rba2=-999. & v_rba2=-999. & w_rba2=-999. & rba2=fltarr(15)-999.
  goto,para_array
endif
fitsg = mpfitpeak(x[bpp],y[bpp], a, nterms = 4, /double, /positive)   ;fitsg=a[3]+a[0]*exp(-0.5*(x-a[1])^2./(a[2]^2))
chisq = (1./(bppc - 4)) * total(((y[bpp] - fitsg[bpp])/err_ave[bpp])^2)
oplot,x,spline(x[bpp],fitsg,x),color=3L
xyouts,x[1],10^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.92),'i:  '+strmid(strtrim(string(a[0]),2),0,6),charsize=1.5,color=3L
xyouts,x[1],10^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.86),'v:  '+strmid(strtrim(string((a[1]-wave0)/wave0*3e5),2),0,4),charsize=1.5,color=3L
xyouts,x[1],10^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.80),'w:  '+strmid(strtrim(string(a[2]/wave0*3e5*sqrt(2)),2),0,4),charsize=1.5,color=3L


;RBp profile and RBp parameters
ymax=max(smooth(interpol(y,nw*100,/spline),3),sub_m)  &  xvec=interpol(x,nw*100,/spline)
x_rel = interpol(x,nw*10,/spline)-xvec[sub_m]
rbstr = gen_rb_profile(x_rel/wave0*2.999e5, interpol(y,nw*10,/spline), steps, 20.)
rba2= (rbstr.red - rbstr.blue)/ymax  &  rba2=smooth(rba2,3)
i_rba2=min(interpol(rba2,n_elements(steps)*10,/spline),sub_m)  
v_rba2=(sub_m/10.+1)*10 ;use the minimum of the R-B profile to determine the velocity
sub_e=where(interpol(rba2,n_elements(steps)*10,/spline) le i_rba2/2.72,sub_n) 
if sub_e[0] eq -1 then w_rba2=-999. else w_rba2=(sub_e[sub_n-1]-sub_e[0])/2  ;exponential width
plot,steps+10,rba2,title='RBp profile',position=posi+[0.6,0,0.38,0],/noerase,$
    yrange=[-0.5,0.5],ystyle=1,xtitle=xtitleRB,ytitle='RB Asymmetry (percentage)',charsize=1.5,/nodata
oplot,steps+10,rba2,color=2L
xyouts,steps[1],(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.92),'i:  '+strmid(strtrim(string(i_rba2),2),0,4),charsize=1.5,color=2L
xyouts,steps[1],(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.86),'v:  '+strmid(strtrim(string(v_rba2),2),0,4),charsize=1.5,color=2L
xyouts,steps[1],(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.80),'w:  '+strmid(strtrim(string(w_rba2),2),0,4),charsize=1.5,color=2L

para_array:
para_sgf[0,n,m]=a[0] ;peak intensity
para_sgf[1,n,m]=(a[1]-wave0)/wave0*3e5 ;Doppler shift
para_sgf[2,n,m]=a[2]/wave0*3e5*sqrt(2) ;exponential width
para_sgf[3,n,m]=a[3] ;background
chisq_sgf[n,m]=chisq
para_rb[0,n,m]=i_rba2
para_rb[1,n,m]=v_rba2
para_rb[2,n,m]=w_rba2
rb_array[n,m,*]=rba2

fit_para[0,n,m]=a[0]
fit_para[1,n,m]=a[1]
fit_para[2,n,m]=a[2]

if chisq ne -999. then begin
r_array[n,m,*]=rbstr.red/ymax
b_array[n,m,*]=rbstr.blue/ymax
endif else begin
r_array[n,m,*]=0
b_array[n,m,*]=0
endelse
wait,1.e-6

endelse
endfor
print,'y=',m,';ny=',ny

endfor

filesave= str_replace(file_lps,'datacube','sgf_rbp_new') 
save,filename=filesave,para_sgf,chisq_sgf,para_rb,rb_array,r_array,b_array,lineid,date_obs,solar_x,solar_y,wave0,time,fit_para
endfor

end
