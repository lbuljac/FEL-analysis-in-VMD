package provide global_fel_analysis 1.0

namespace eval ::GlobalFEL:: {
    variable selection_text "protein and name CA"
    variable eps 0.05
    variable min_pts 53
    variable temperature 300
    variable ref_id 0
}

proc ::GlobalFEL::calculate_metrics {mol_id ref_mol_id sel_str output_file} {
    set sel ""
    set ref ""
    set outfile ""
    
    if {[catch {
        set sel [atomselect $mol_id $sel_str]
        set ref [atomselect $ref_mol_id $sel_str frame 0]
        set nf [molinfo $mol_id get numframes]

        set outfile [open $output_file "w"]
        puts $outfile "Frame RMSD Rg"

        for {set i 0} {$i < $nf} {incr i} {
            $sel frame $i
            set trans_mat [measure fit $sel $ref]
            $sel move $trans_mat

            set rmsd [measure rmsd $sel $ref]
            set rg [measure rgyr $sel]

            puts $outfile "$i $rmsd $rg"
        }
    } err]} {
        if {$outfile ne ""} { catch {close $outfile} }
        if {$sel ne ""} { catch {$sel delete} }
        if {$ref ne ""} { catch {$ref delete} }
        error "Error while calculating for molecule $mol_id: $err"
    }

    if {$outfile ne ""} { catch {close $outfile} }
    if {$sel ne ""} { catch {$sel delete} }
    if {$ref ne ""} { catch {$ref delete} }
}

proc ::GlobalFEL::ensure_metrics_calculated {} {
    set ref_id $::GlobalFEL::ref_id
    
    if {![file exists "md1_metrics.txt"] || ![file exists "md2_metrics.txt"] || ![file exists "md3_metrics.txt"]} {
        puts "(VMD Plugin) Calculating RMSD and Rg data for 3 MD simulation replicas..."
        ::GlobalFEL::calculate_metrics 1 $ref_id $::GlobalFEL::selection_text "md1_metrics.txt"
        ::GlobalFEL::calculate_metrics 2 $ref_id $::GlobalFEL::selection_text "md2_metrics.txt"
        ::GlobalFEL::calculate_metrics 3 $ref_id $::GlobalFEL::selection_text "md3_metrics.txt"
    } else {
        puts "(VMD Plugin) Pre-calculated RMSD/Rg data are used."
    }
}

proc ::GlobalFEL::run_knn_only {} {
    if {[catch { ::GlobalFEL::ensure_metrics_calculated } err]} {
        puts "Error while preparing metrics: $err"
        return
    }

    set script_path [file join [file dirname [info script]] "scripts" "fel_analysis.py"]
    puts "(VMD Plugin) Generating k-NN distance plot with eps = $::GlobalFEL::eps i minPts = $::GlobalFEL::min_pts..."
    
    if {[catch {
        exec cmd /c python $script_path "knn" "md1_metrics.txt" "md2_metrics.txt" "md3_metrics.txt" $::GlobalFEL::eps $::GlobalFEL::min_pts $::GlobalFEL::temperature
    } err]} {
        puts "Error while generating k-NN graph: $err"
        return
    }

    puts "(VMD Plugin) k-NN graph saved in 'kNN_distance_plot.png'."
}

proc ::GlobalFEL::run_full_pipeline {} {
    if {[catch { ::GlobalFEL::ensure_metrics_calculated } err]} {
        puts "Error while preparing metrics: $err"
        return
    }

    set script_path [file join [file dirname [info script]] "scripts" "fel_analysis.py"]
    puts "(VMD Plugin) Starting full FEL analysis..."
    
    if {[catch {
        exec cmd /c python $script_path "full" "md1_metrics.txt" "md2_metrics.txt" "md3_metrics.txt" $::GlobalFEL::eps $::GlobalFEL::min_pts $::GlobalFEL::temperature
    } err]} {
        puts "Error during full FEL analysis: $err"
        return
    }

    puts "(VMD Plugin) Full FEL analysis completed! All plots and statistics are generated successfully."
}

proc ::GlobalFEL::create_gui {} {
    set w .global_fel_window

    if {[winfo exists $w]} {
        wm deiconify $w
        raise $w
        focus -force $w
        return
    }

    toplevel $w
    wm title $w "Global FEL Analysis"

    # Atom selection
    frame $w.f_sel
    label $w.f_sel.l -text "Atom Selection:"
    entry $w.f_sel.e -textvariable ::GlobalFEL::selection_text -width 35
    pack $w.f_sel.l $w.f_sel.e -side left -padx 5 -pady 5
    pack $w.f_sel -fill x

    # DBSCAN parameters
    frame $w.f_dbscan
    label $w.f_dbscan.l_eps -text "eps:"
    entry $w.f_dbscan.e_eps -textvariable ::GlobalFEL::eps -width 6
    label $w.f_dbscan.l_pts -text "minPts:"
    entry $w.f_dbscan.e_pts -textvariable ::GlobalFEL::min_pts -width 6
    pack $w.f_dbscan.l_eps $w.f_dbscan.e_eps $w.f_dbscan.l_pts $w.f_dbscan.e_pts -side left -padx 5
    pack $w.f_dbscan -fill x -pady 5

    # Setting the GUI button
    frame $w.f_btn
    button $w.f_btn.btn_knn -text "Set DBSCAN parameters" -command ::GlobalFEL::run_knn_only
    button $w.f_btn.btn_run -text "Run Analysis" -command ::GlobalFEL::run_full_pipeline
    pack $w.f_btn.btn_knn $w.f_btn.btn_run -side left -padx 5 -pady 10
    pack $w.f_btn -anchor center

    raise $w
    focus -force $w
}

# Registration to Extensions -> Analysis GUI menu
if {[info exists ::vmd_version]} {
    vmd_install_extension global_fel_analysis ::GlobalFEL::create_gui "Analysis/Global FEL Analysis"
}