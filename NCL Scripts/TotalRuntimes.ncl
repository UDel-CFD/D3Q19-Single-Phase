; NCL Code for plotting total run times with linear regression

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;Declare processor count, and timestep counts
    numprocesses = 400
    runtimes = (/50, 100, 1000, 2000/)
    timedata = new((/numprocesses,dimsizes(runtimes)/), "float")
    plotdata = new((/2, max(runtimes)/), "float")

;File directory
    fileroot = systemfunc("pwd") + "/Data/runtimes"
    filenames = fileroot + runtimes(:) + ".dat"

; Read in runtime files
    do i=0, dimsizes(runtimes)-1
        filedata = asciiread(filenames(i), (/numprocesses,2/), "float") ;Remember ncl is [row, col]
        timedata(:,i) = filedata(:,1)
    end do

; Average/stddev the time data for all processors
    averagedata = dim_avg_n(timedata,0)
    stddevdata = dim_stddev_n(timedata,0)

    print(sprintf("%12.9f",averagedata))
    print(sprintf("%12.9f",stddevdata*2))

    plotdata(0,runtimes(:)-1) = averagedata(:)

; Linear Regression
; Use regline_stats for more data on the fit
    xspan =ispan(0,max(runtimes)-1,1)
    rc = regline_stats(int2flt(runtimes), averagedata) ;Arrays must be of same type!
    plotdata(1,:) = rc*(xspan) + rc@yintercept ;mx + b
    print(sprintf("%12.9f",rc))
; =================================
; Using gsun so we don't have to worry about app, workstations, data, etc.
    wks   = gsn_open_wks ("x11","xy")

; Set up plot format and plot data sets
    res                   = True                                    ; plot mods desired
    res@gsnMaximize       = True                                    ; maximize plot
    res@tiMainString      = "Main Loop Wall-clock Time"             ; add title
    res@gsnCenterString   = "~Z70~nx=199, ny=400, nz=200, 400 MPI Tasks, "\
                                +"npart=270"                        ; add center label for description, ~z#~ changes font
    res@tiMainJust        ="CenterCenter"                           ; Set title position
    res@tiMainFontHeightF = 0.025                                   ; set title font
    res@tiXAxisString     = "Time Steps"                            ; X Axis label
    res@tiYAxisString     = "Wall-clock time (s)"                   ; Y axis label

    res@tmXBMode        = "Explicit"                                ; specify x-axis tick labels
    res@tmXBValues      = ispan(0,max(runtimes)+100,200)            ; set x major tick spacing
    res@tmXBLabels      = "" + res@tmXBValues                       ; set x major tick labels
    res@tmXBMinorValues = ispan(0,max(runtimes)+100,50)             ; set x minor tick spacing
    res@tmXBLabelFontHeightF = 0.013                                ; set font size for x labels, default = 0.02
    res@trYMinF               = 0
    res@trYMaxF               = max(plotdata(0,:))+25

    res@xyDashPattern     = 0                                       ; choose solid line
    res@xyMarkLineModes   = (/"Markers","Lines"/)                   ; choose to use markers or line
    res@xyMarkers         =  16                                     ; choose type of marker  
    res@xyMarkerColor     = "cyan"                                  ; Marker color
    res@xyLineColor       = "red"                                   ; line color
    res@xyMarkerSizeF     = 0.01                                    ; Marker size (default 0.01)
    res@xyLineThicknessF  = 0.5                                     ; line thickness

    res@gsnFrame        = False                                     ; don't advance frame yet
    res@gsnDraw         = False                                     ; don't draw plot

    plot  = gsn_csm_xy (wks,xspan,plotdata,res)                     ; create plot

; =================================
; Add linear regression label
    txres               = True                     ; text mods desired
    txres@txFontHeightF = 0.018                    ; font smaller. default big

    text = "y = " + sprintf("%12.9f",rc) + "x + " + sprintf("%12.9f",rc@yintercept) +\
                "~C~r = "+sprintf("%12.9f",rc@r2)
    txt_summer = gsn_add_text(wks,plot,text, 1200, 50.0, txres)

; =================================
; Add error bars 2*stddev for 95%
    polyres                   = True                       ; marker resources
    polyres@gsMarkerIndex     = 1                          ; index for circle
    polyres@gsMarkerSizeF     = .02                        ; size
    polyres@gsLineThicknessF  = 3                          ; line thickness

    error_bar = new(dimsizes(runtimes),graphic)
    centers   = new(dimsizes(runtimes),graphic) 
; loop through the points
    do t=0,dimsizes(runtimes)-1
; add a polymarker as the center of the error bar located at each data point
        centers(t)   = gsn_add_polymarker(wks,plot,runtimes(t),averagedata(t),polyres)
; the error bar is draw in two segements, above and below the true value of u(t)
        err = stddevdata(t)*2
        error_bar(t) = gsn_add_polyline(wks,plot,(/runtimes(t),runtimes(t)/),(/averagedata(t)+\
                                        err,averagedata(t)-err/),polyres)
    end do

; =================================
; Draw plot!
    draw(plot)
    frame(wks)
end