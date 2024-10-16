function subset_selected = downsample(target_sample_size,vpop_inclusion_bool)

    %Downsample a virtual population (i.e. selected via Allen/Reiger alg.)
    %into a specific sample size that must be smaller than the original
    %population. Uses equi-probability resampling. Assumes the virtual
    %population itself is specified as a boolean vector that is the same
    %length as the plausible population from which it was selected. This
    %function will select a subset of the true/non-zero elements from this
    %passed Vpop subset as specified in the boolean vector.

    %Input
    %target_sample_size -- interger specifying target sample size, must be
    %                       smaller than sum(vpop_inclusion_bool) (i.e.
    %                       size of vpop)
    %vpop_inclusion_bool -- boolean vector, same length as original
    %                       plausible population from which Vpop was
    %                       selected, 1 indicates plausible was in Vpop, 0
    %                       indicates plausible was not included. This
    %                       algorithm will select a subset of indices
    %                       corresponding to non-zero values.

    %Output
    %subset_selected -- boolean, same size as original plausible population
    %                   with 1's indicated inclusion in final downsampled
    %                   Vpop. These non-zero values will always be a subset
    %                   of vpop_inclusion_bool non-zeros.

    %get number of plausible patients (default to 10k for paper)
    num_plausibles = length(vpop_inclusion_bool);
    %sum to get number of virtual patients selected by Allen/Riger Vpop
    num_virtuals = sum(vpop_inclusion_bool);
    %list of plausible patient indices
    plausible_ids = 1:num_plausibles;
    %subset plausible ids based on inclusion in Allen/Rieger Vpop
    plausible_ids_selected = plausible_ids(vpop_inclusion_bool==1);
    
    %seed random number generator for reproducibility
    %rng(1,'twister') 
    %select a subset of virtual patient ids from Allen/Riger selection to
    %downsample for target sample size
    virtual_subset_indices = randsample(num_virtuals,target_sample_size);
    
    %subset the selected plausibles ids to get downsampled subset
    plausible_ids_virtual_subset = plausible_ids_selected(virtual_subset_indices);
    
    %create a new selected vector (same length as plausibles) that selects
    %downsampled subset of virtual patients orginally selected by Allen/Rieger
    subset_selected = zeros(num_plausibles,1);
    subset_selected(plausible_ids_virtual_subset) = 1;
    subset_selected = subset_selected == 1;
end