function weighted_sample = subsamp(sample_size,weight)

    %A helper function to sample weith replacment specific number of 
    %plausible patients according to probability given in weight vector. 

    %Input
    % sample_size -- desired sample size
    % weight -- weight vector, same length as number of plausible patients,
    %           where weight is proportional to desired probability of
    %           selection

    %sample w/ replacment indices from 1 to length weight according to 
    % probability in weight vector
    rs=randsample(length(weight),sample_size,true,weight);
    
    %create vector same length as weight vector to store resampled counts
    weighted_sample=zeros(length(weight),1);
    %count the number of times each index was sampled, store in vector
    for i=1:length(rs)
        weighted_sample(rs(i))=weighted_sample(rs(i))+1; 
    end   
end