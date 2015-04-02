; NCL Code for plotting subroutine runtimes

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;Declare processor, subroutine and timestep counts
    numprocesses = 400
    numtimesteps = 2000
    numsubroutines = 9
    subroutinelabels = (/"MRT Collison", "Beads Collision", "Streaming", "Beads Lubforce"\
                        ,"Beads Move", "Beads Redistribute", "Beads Links", "Beads Filling"\
                        ,"Macrovar"/)

    processes = ispan(0, numprocesses-1, 1)
    timedata = new((/numprocesses, numtimesteps, numsubroutines/), "float")

;File directory
    fileroot = systemfunc("pwd") + "/Data/Subroutines/bmmatlab."
    filenames = fileroot + sprinti("%0.3i",processes(:))

; Read in subroutine files
    do i=0, numprocesses-1
        filedata = asciiread(filenames(i), (/numtimesteps,numsubroutines/), "float") ;Remember ncl is [row, col]
        timedata(i,:,:) = filedata(:,:)
    end do

; Averaging and stddev times
    timeavgtemp = dim_avg_n(timedata,1) ; Average between timesteps
    averagedata = dim_avg_n(timeavgtemp, 0) ; Average between processes
    stddevdata = dim_stddev_n(timeavgtemp, 0) ; Stddev between processes

    print(sprintf("%12.9f",averagedata))
    print(sprintf("%12.9f",stddevdata*2))

    x = ispan(1, numsubroutines, 1)
;=========================================
;Open gsn workstation
    wks = gsn_open_wks("x11","bar")

; Set up plot format and data sets
    res                       = True                                    ; plot mods desired
    res@tiMainString          = "Subroutine Wall-Clock Times"           ; plot title
    res@tiYAxisString         = "Wall-clock time (s)"                   ; Y axis label 
    res@gsnCenterString   = "~Z70~nx=199, ny=400, nz=200, 400 MPI Tasks, "\
                                +"npart=270"                            ; add center label for description, ~z#~ changes font

    res@gsnXYBarChart         = True                                    ; turn on bar chart
    res@gsnXYBarChartBarWidth = 0.75                                    ; change bar widths
    res@gsnXYBarChartColors = (/"firebrick","red","orange","green", \
    "navy","blue","SkyBlue","SlateBlue", "orchid3"/)                    ; choose colors
  
    res@tmXBOn                = False                                   ; turn off tickmarks at bot
    res@trYMinF               = 0                                       ; bring bars down to zero
    res@trXMinF               = 0                                       ; adds space on either end
    res@trXMaxF               = 10                                      ; of the 1st and last bars

    res@gsnFrame        = False                                         ; don't advance frame yet
    res@gsnDraw         = False                                         ; don't draw plot

    plot = gsn_csm_xy (wks,x,averagedata,res)                           ; create plot

; =======================================
; Add bar labels
    txres                     = True                 ; text mods desired
    txres@txFontHeightF       = 0.018                ; default size is HUGE!
    txres@txAngleF            = 52.                  ; text angle
    txres@txJust              = "TopRight"           ; moves text anchor to top right
  
    do n = 1, numsubroutines
        gsn_text(wks,plot,subroutinelabels(n-1),n-0.018,0,txres) ; add labels
    end do

; =================================
; Add error bars 2*stddev for 95%
    polyres                   = True                       ; marker resources
    polyres@gsMarkerIndex     = 1                          ; index for circle
    polyres@gsMarkerSizeF     = .02                        ; size
    polyres@gsLineThicknessF  = 3                          ; line thickness

    error_bar = new(numsubroutines,graphic)
    centers   = new(numsubroutines,graphic) 
; loop through bars
    do t=0,numsubroutines-1
; add a polymarker as the center of the error bar located at each data point
        centers(t)   = gsn_add_polymarker(wks,plot,t+1,averagedata(t),polyres)
; the error bar is draw in two segements, above and below the true value of u(t)
        err = stddevdata(t)*2
        error_bar(t) = gsn_add_polyline(wks,plot,(/t+1,t+1/),(/averagedata(t)+\
                                        err,averagedata(t)-err/),polyres)
    end do

; =======================================
; Draw plot!
    draw(plot)
    frame(wks)
end
