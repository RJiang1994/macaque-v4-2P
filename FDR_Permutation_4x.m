% The mean df/f0 map and the variance map can be found in Data_4x_A.mat and Data_4x_B.mat
% The df/f0 map for each single trial can be found in osf.io/xmj5c/


rng('shuffle')

load DeltaFluoreFilteredB.mat
ClusterSizeCV = zeros(1,10000);
ClusterSizeCN = zeros(1,10000);
repeat = 18; %10

tic;

for k = 1:10000
        k
        
        % for fdr map
        CV = 1:32*repeat;
        CN = 32*repeat+1:(32+48)*repeat;
        BAR = (32+48)*repeat+1:(32+48+16)*repeat;
              
        % random permutation
%         relabel = (1:length(96*repeat));  % for fdr map
        relabel = randperm(96*repeat);  % for permutation
        cv = relabel(CV);
        cn = relabel(CN);
        bar = relabel(BAR);
        
        % dff = df/f0 smoothed
        dff_cv = dff(:,:,cv);
        dff_cn = dff(:,:,cn);
        dff_bar = dff(:,:,bar);
        
        % average response map
        mdff_cv = mean(dff_cv,3);
        mdff_cn = mean(dff_cn,3);
        mdff_bar = mean(dff_bar,3);
       
        var_cv = var(dff_cv,0,3);
        var_cn = var(dff_cn,0,3);
        var_bar = var(dff_bar,0,3);
     
        % group variance map
        var_cvcn = ( (length(cv)-1)*var_cv + (length(cn)-1)*var_cn ) / ( length(cv) + length(cn) - 2 );
        var_cvbar = ( (length(cv)-1)*var_cv + (length(bar)-1)*var_bar ) / ( length(cv) + length(bar) - 2 );
        var_cnbar = ( (length(cn)-1)*var_cn + (length(bar)-1)*var_bar ) / ( length(cn) + length(bar) - 2 );
        
        % t-statistic map
        t_cvcn = ( mdff_cv - mdff_cn ) ./ sqrt(var_cvcn*(1/length(cv)+1/length(cn)));
        t_cvbar = ( mdff_cv - mdff_bar ) ./ sqrt(var_cvbar*(1/length(cv)+1/length(bar)));
        t_cnbar = ( mdff_cn - mdff_bar ) ./ sqrt(var_cnbar*(1/length(cn)+1/length(bar)));
              
        p_cvcn = 1 - tcdf(t_cvcn,length(cv)+length(cn)-2);
        p_cvbar = 1 - tcdf(t_cvbar,length(cv)+length(bar)-2);
        p_cnbar = 1 - tcdf(t_cnbar,length(cn)+length(bar)-2);

        % p value map
        p_cv = max(p_cvcn,p_cvbar);
        p_cn = max(1-p_cvcn,p_cnbar);
     
        % fdr map
        fdr_cv = mafdr(p_cv(:),'BHFDR',1);
        
%         h_cv = find(fdr_cv<0.01);
        h_cv = find(p_cv(:)<0.01);
        bw = zeros(512);
        bw(h_cv) = 1;
        figure,imshow(bw);
        bw1=bw;
        CC_cv = bwconncomp(bw1,4);
        sizelist = zeros(1,CC_cv.NumObjects);
        for ci = 1:CC_cv.NumObjects
            sizelist(ci) = length(CC_cv.PixelIdxList{ci});
        end
        if ~isempty(sizelist)    
            ClusterSizeCV(k) = max(sizelist);
        end
        
        % fdr map
        fdr_cn = mafdr(p_cn(:),'BHFDR',1);
        
%         h_cn = find(fdr_cn<0.01);
        h_cn = find(p_cn(:)<0.01);
        bw = zeros(512);
        bw(h_cn) = 1;
        figure,imshow(bw);
        bw2=bw;
        CC_cn = bwconncomp(bw2,4);
        sizelist = zeros(1,CC_cn.NumObjects);
        for ci = 1:CC_cn.NumObjects
            sizelist(ci) = length(CC_cn.PixelIdxList{ci});
        end
        if ~isempty(sizelist)
            ClusterSizeCN(k) = max(sizelist);
        end
  
end

toc;
