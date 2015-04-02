; NCL Code for plotting subroutine runtimes

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;Declare processor, subroutine and timestep counts
    numprocesses = 400
    numtimesteps = 2000
    numsubroutines = 9
    subroutine = 5 ; Which sub routine you want to compare, see labels below

    subroutinelabels = (/"MRT Collison", "Beads Collision", "Streaming", "Beads Lubforce"\
                        ,"Beads Move", "Beads Redistribute", "Beads Links", "Beads Filling"\
                        ,"Macrovar"/)

    processes = ispan(0, numprocesses-1, 1)
    timedata = new((/numprocesses, numtimesteps, numsubroutines/), "float")
    timedata0 = new((/numprocesses, numtimesteps, numsubroutines/), "float") ; holdsold time data
    plotdata = new((/4, numprocesses/), "float") ; 4 sets of data, 2 scatter, 2 average lines

;New File directory
    fileroot = systemfunc("pwd") + "/Data/Subroutines/bmmatlab."
    filenames = fileroot + sprinti("%0.3i",processes(:)) ; Pad with 0s

;Old File directory
    fileroot = systemfunc("pwd") + "/Data/Subroutines/old/bmmatlab."
    filenames0 = fileroot + sprinti("%0.3i",processes(:)) ; Pad with 0s

; Read in subroutine files
    do i=0, numprocesses-1
        filedata = asciiread(filenames(i), (/numtimesteps,numsubroutines/), "float") ;Remember ncl is [row, col]
        timedata(i,:,:) = filedata(:,:)
    end do

    do i=0, numprocesses-1
        filedata = asciiread(filenames0(i), (/numtimesteps,numsubroutines/), "float") ;Remember ncl is [row, col]
        timedata0(i,:,:) = filedata(:,:)
    end do


; Averaging wall times
    averagedata = dim_avg_n(timedata,1) ; Average between timesteps
    averagedata0 = dim_avg_n(timedata0, 1) ; Average between timesteps
    plotdata(0,:) = averagedata(:,subroutine)
    plotdata(2,:) = averagedata0(:,subroutine)

    processavg = dim_avg_n(averagedata, 0)
    processavg0 = dim_avg_n(averagedata0, 0)
    plotdata(1,:) = processavg(subroutine) ; Average between processes
    plotdata(3,:) = processavg0(subroutine) ; Average between processes

    print("Original:")
    print(sprintf("%12.9f",processavg0(subroutine)))
    print("Optimized:")
    print(sprintf("%12.9f",processavg(subroutine)))
;=========================================
;Open gsn workstation
    wks = gsn_open_wks("x11","xy")

; Set up plot format and data sets
    res                       = True                                            ; plot mods desired
    res@tiMainString          = subroutinelabels(subroutine)+\
                                    " Wall-Clock Times"                         ; plot title
    res@tiXAxisString         = "MPI Task Number"                            ; Y axis label                                 
    res@tiYAxisString         = "Wall-clock time (s)"                           ; X axis label
    res@gsnCenterString   = "~Z70~nx=199, ny=400, nz=200, 400 MPI Tasks, "\
                                +"npart=270"                                    ; add center label for description, ~z#~ changes font

    res@xyDashPattern     = 0                                                   ; choose solid line
    res@xyMarkLineModes   = (/"Markers","Lines","Markers","Lines"/)             ; choose to use markers or line
    res@xyMarkers         = (/0,0/)                                             ; choose type of marker  
    res@xyMarkerColors    = (/"Red","Red","Blue","Blue"/)                       ; Marker color
    res@xyLineColors      = (/"Red","Red","Blue","Blue"/)                       ; Line colors, Must assing markers line colors too
    res@xyMarkerSizes     = (/0.01,0.01/)                                       ; Marker size (default 0.01)
    res@xyLineThicknessF  = 2                                                   ; line thickness

    res@pmLegendDisplayMode    = "Always"                                       ; turn on legend
    res@pmLegendSide           = "Right"                                        ; Change location of
    res@gsnMaximize            = True                                           ; Make window include all of legend
    res@pmLegendParallelPosF   = 0.66                                            ; move units right
    res@pmLegendOrthogonalPosF = 0.1                                              ; move units down 
    res@pmLegendWidthF         = 0.15                                           ; Change width and
    res@pmLegendHeightF        = 0.18                                           ; height of legend.
    res@lgLabelFontHeightF     = .02                                            ; label font height
    res@xyExplicitLegendLabels = (/"Optimized","Optimized Average",\
                                  "Original","Orginal Average"/)                ; create explicit labels

    res@gsnFrame        = False                                                 ; don't advance frame yet
    res@gsnDraw         = False                                                 ; don't draw plot

    plot = gsn_csm_xy (wks,processes,plotdata,res)                              ; create plot

; =======================================
; Draw plot!
    draw(plot)
    frame(wks)
end
