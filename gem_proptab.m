function out_table = gem_proptab(model, blockedRxns)
% Function to calculate a table of GEM properties 
% https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialModelProperties.html
% Input
%   struc model:    A cobra model
%   logic blockedRxns:  if true time consuming FVA is run to find blocked
%                       RXNs (default:false)
if ~exist('blockedRxns', 'var')
    blockedRxns=false;
end
clear TableProp
r=1;
TableProp(r,:) = {'Model'}; r=r+1;

TableProp(r, 1) = {'Reactions'};
%Number of reactions
TableProp{r, 2} = num2str(length(model.rxns));
r = r + 1;

TableProp(r, 1) = {'Metabolites'};
%Number of metabolites
TableProp{r, 2} = num2str(length(model.mets));
r = r + 1;

TableProp(r, 1) = {'Metabolites (unique)'};
%Metabolites contain compartment in Name enclosed by brackets
%Code removes it and counts unique names after that
[g, remR3M] = strtok(model.mets,'[');
TableProp{r, 2} = num2str(length(unique(g)));
r = r + 1;
%Following the same annotation unique compartments are extracted
TableProp(r, 1) = {'Compartments (unique)'};
TableProp{r, 2} = num2str(length(unique(remR3M)));
r = r + 1;

TableProp(r, 1) = {'Deadends'};
D3M = detectDeadEnds(model);
TableProp{r, 2} = num2str(length(D3M));
r = r + 1;

TableProp(r, 1) = {'Size of S'};
TableProp{r, 2} = strcat(num2str(size(model.S,1)),'; ',num2str(size(model.S,2)));
r = r + 1;
 %rank (linearly independent columns/rows)
TableProp(r, 1) = {'Rank of S'};
TableProp{r, 2} = strcat(num2str(rank(full(model.S))));
r = r + 1;
%nonzero elements in S
TableProp(r, 1) = {'Percentage nz'};
TableProp{r, 2} = strcat(num2str((nnz(model.S)/(size(model.S,1)*size(model.S,2)))));
r = r + 1;

%GPR rules 
TableProp(r,1)= {'Assoc. GPR'};
TableProp{r,2}= num2str(sum(~(cellfun(@isempty,model.rules))));
r=r+1;

%unique genes
TableProp(r,1)={'Genes (unique)'};
TableProp{r,2}=num2str(length(unique(model.genes)));
r=r+1;
%blocked reactions
if blockedRxns
    TableProp(r,1)={'Blocked Rxns'};
    try
        TableProp{r,2}=num2str(length(findBlockedReaction(model)));
    catch
        warning('FVA for blocked Rxns not possible\n setting Blocked Rxns to NaN')
        TableProp{r,2}='NaN';
    end
end

%view table
out_table=TableProp;
end
