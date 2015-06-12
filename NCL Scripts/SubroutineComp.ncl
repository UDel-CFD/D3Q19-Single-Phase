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
    timedata = new((/numprocesses, numtimesteps, numsubroutines,2/), "float")

;File directories
    fileroot = systemfunc("pwd") + "/Data/Subroutines/bmmatlab."
    filenamesnew = fileroot + sprinti("%0.3i",processes(:))
    fileroot = systemfunc("pwd") + "/Data/Subroutines/Old/bmmatlab."
    filenamesold = fileroot + sprinti("%0.3i",processes(:))

; Read in subroutine files
    do i=0, numprocesses-1
        filedata = asciiread(filenamesnew(i), (/numtimesteps,numsubroutines/), "float") ;Remember ncl is [row, col]
        timedata(i,:,:,0) = filedata(:,:)
    end do
    do i=0, numprocesses-1
        filedata = asciiread(filenamesold(i), (/numtimesteps,numsubroutines/), "float") ;Remember ncl is [row, col]
        timedata(i,:,:,1) = filedata(:,:)
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
    res@gsnXYBarChartBarWidth = 0.40                                    ; change bar widths
    res@tmXBValues            = ispan(1,numsubroutines,1)
  
    res@tmXBOn                = False                                   ; turn off tickmarks at bot
    res@trYMinF               = 0                                       ; bring bars down to zero
    res@trXMinF               = 0                                       ; adds space on either end
    res@trXMaxF               = 10                                      ; of the 1st and last bars

    res@gsnFrame        = False                                         ; don't advance frame yet
    res@gsnDraw         = False                                         ; don't draw plot

    res@gsnXYBarChartColors = (/"lightblue"/)
    plot1 = gsn_csm_xy (wks,fspan(.80,numsubroutines-.20,numsubroutines),averagedata(:,0),res)  ; create plot
    res@gsnXYBarChartColors = (/"blue"/)
    plot2 = gsn_csm_xy (wks,fspan(1.2,numsubroutines+.20,numsubroutines),averagedata(:,1),res)  ; create plot

; =======================================
; Add bar labels
    txres                     = True                 ; text mods desired
    txres@txFontHeightF       = 0.018                ; default size is HUGE!
    txres@txAngleF            = 52.                  ; text angle
    txres@txJust              = "TopRight"           ; moves text anchor to top right
  
    do n = 1, numsubroutines
        gsn_text(wks,plot1,subroutinelabels(n-1),n-0.018,0,txres) ; add labels
    end do
; =======================================
; Add color labels
    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.3           ; labelbar width
    lbres@vpHeightF          = 0.1           ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
    lbres@vpXF               = 0.50          ; labelbar position in x direction 
    lbres@vpYF               = 0.80          ; labelbar position in y direction
    lbres@lbFillColors       = (/"lightblue","blue"/)
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = 0.018         ; font height
    lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    lbres@lbPerimOn          = False
    lbres@lgPerimColor      = "white"
    labels = (/"Optimized","Original"/)
    lbid =  gsn_create_labelbar(wks,2,labels,lbres) ; draw right labelbar column
; =======================================
; Draw plot!
    draw(plot1)
    draw(plot2)
    draw(lbid)
    frame(wks)
end
