addpath(genpath('/Users/christopher_penfold/Desktop/Code/SCUBA_/drtoolbox/'))

RNAseq_preprocess('PGC2',1)
SCUBA('PGC2')

D1 = importdata('/Users/christopher_penfold/Desktop/BranchingGPs/results/primordial_germ_cells/OtherPseudotime/SCUBA/PGC2Data.txt')

ps = load('/intermediate_files/PGC2PData.mat');
cells = D1.textdata(1,2:end);
pseudotime = ps.pro.pseudotime;

fid = fopen('./Pseudotime.csv','a')
for i = 1:length(pseudotime)
    fprintf(fid,'%s,',cells{i})
    fprintf(fid,'%f\n',pseudotime(i))
end
fclose(fid)