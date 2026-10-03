function tf = prereq_satisfied(net, i, connected)
p = net.parents{i};
if isempty(p)
    tf = true;
elseif net.parentLogic(i) == "OR"
    tf = any(connected(p));
else
    tf = all(connected(p));
end
end
