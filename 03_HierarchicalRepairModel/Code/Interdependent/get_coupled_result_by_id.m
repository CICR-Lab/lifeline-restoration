function [scenario,out]=get_coupled_result_by_id(catalog,results,id)
idx=find(strcmp({catalog.id},id),1);
if isempty(idx), error('Coupled scenario ID not found: %s',id); end
scenario=catalog(idx); out=results{idx};
end
