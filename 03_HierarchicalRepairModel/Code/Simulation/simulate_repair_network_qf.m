function [trajectory, diagnostics] = simulate_repair_network_qf(net,nCrews,qF,cfg)

N = net.N;
if ~isscalar(qF) || ~isfinite(qF) || qF < 0 || qF > 1
    error('RepairTask:InvalidFrontierProbability', ...
        'qF must be a finite scalar on [0,1].');
end
nCrews = min(N, max(1, round(nCrews)));
started = false(N,1);
completed = false(N,1);
connected = false(N,1);
repairFinishTime = nan(N,1);
ready = false(N,1);
ready(net.layer==0) = true;
remainingPrerequisites = double(net.parentCount);

priority = net.randomPriority;
[~,globalOrder] = sort(priority,'ascend');
globalPointer = 1;

readyHeap = heap_initialize(N);
initialReady = find(ready);
for j = 1:numel(initialReady)
    task = initialReady(j);
    readyHeap = heap_push(readyHeap, priority(task), task);
end
activeHeap = heap_initialize(N);

currentTime = 0;
nStarted = 0;
nConnected = 0;
nFrontierStarts = 0;
restored = 0;

trajectoryTime = zeros(N+1,1);
trajectoryRestored = zeros(N+1,1);
nTrajectory = 1;
completedBatch = zeros(N,1);
queue = zeros(N,1);

while nConnected < N
    while activeHeap.size < nCrews && nStarted < N
        [readyHeap, frontierTask] = peek_valid_ready( ...
            readyHeap,started);
        while globalPointer <= N && started(globalOrder(globalPointer))
            globalPointer = globalPointer + 1;
        end

        chooseFrontier = frontierTask > 0 && ...
            net.hybridChoice(nStarted+1) < qF;
        if chooseFrontier
            [~, task, readyHeap] = heap_pop(readyHeap);
        else
            task = globalOrder(globalPointer);
            globalPointer = globalPointer + 1;
        end
        started(task) = true;
        nStarted = nStarted + 1;
        nFrontierStarts = nFrontierStarts + double(ready(task));
        finishTime = currentTime + net.duration(task);
        activeHeap = heap_push(activeHeap, finishTime, task);
    end

    if activeHeap.size == 0
        error('RepairTask:SimulationDeadlock', ...
            'No active task remains before all service has reconnected.');
    end

    nextTime = activeHeap.key(1);
    tolerance = cfg.simultaneousTolerance .* max(1,abs(nextTime));
    nBatch = 0;
    while activeHeap.size > 0 && activeHeap.key(1) <= nextTime + tolerance
        [~, task, activeHeap] = heap_pop(activeHeap);
        completed(task) = true;
        repairFinishTime(task) = nextTime;
        nBatch = nBatch + 1;
        completedBatch(nBatch) = task;
    end
    currentTime = nextTime;

    qHead = 1;
    qTail = 0;
    for j = 1:nBatch
        task = completedBatch(j);
        if ready(task) && ~connected(task)
            connected(task) = true;
            nConnected = nConnected + 1;
            restored = restored + net.weight(task);
            qTail = qTail + 1;
            queue(qTail) = task;
        end
    end

    while qHead <= qTail
        parent = queue(qHead);
        qHead = qHead + 1;
        childList = net.children{parent};
        for h = 1:numel(childList)
            child = childList(h);
            if ready(child)
                continue;
            end
            switch net.dependencyLogic
                case {'single','and2'}
                    remainingPrerequisites(child) = ...
                        remainingPrerequisites(child) - 1;
                    becomesReady = remainingPrerequisites(child) == 0;
                case 'or2'
                    becomesReady = true;
                otherwise
                    error('RepairTask:UnknownDependencyLogic', ...
                        'Unknown dependency logic: %s', net.dependencyLogic);
            end
            if becomesReady
                ready(child) = true;
                if completed(child) && ~connected(child)
                    connected(child) = true;
                    nConnected = nConnected + 1;
                    restored = restored + net.weight(child);
                    qTail = qTail + 1;
                    queue(qTail) = child;
                elseif ~started(child)
                    readyHeap = heap_push(readyHeap, priority(child), child);
                end
            end
        end
    end

    if restored > trajectoryRestored(nTrajectory) || ...
            (nConnected==N && currentTime>trajectoryTime(nTrajectory))
        sameTime = currentTime <= trajectoryTime(nTrajectory) + ...
            cfg.simultaneousTolerance.*max(1,abs(currentTime));
        if sameTime
            trajectoryRestored(nTrajectory) = restored;
        else
            nTrajectory = nTrajectory + 1;
            trajectoryTime(nTrajectory) = currentTime;
            trajectoryRestored(nTrajectory) = restored;
        end
    end
end

trajectory.time = trajectoryTime(1:nTrajectory);
trajectory.restored = trajectoryRestored(1:nTrajectory);

diagnostics.chiFrontier = nFrontierStarts ./ N;
diagnostics.crewUtilization = sum(net.duration) ./ (nCrews .* currentTime);
diagnostics.maxRestorationJump = max(diff(trajectory.restored));
diagnostics.finalService = restored;
diagnostics.isMonotone = all(diff(trajectory.restored) >= -cfg.numericTolerance);
diagnostics.allReconnected = all(connected) && nConnected==N;
diagnostics.repairFinishTime = repairFinishTime;

if ~diagnostics.allReconnected || ~diagnostics.isMonotone || ...
        abs(diagnostics.finalService-1) > 100*cfg.numericTolerance
    error('RepairTask:InvalidTrajectory', ...
        'Generated restoration trajectory failed validity checks.');
end
end

function [heap, task] = peek_valid_ready(heap,started)
task = 0;
while heap.size > 0
    candidate = heap.task(1);
    if ~started(candidate)
        task = candidate;
        return;
    end
    [~,~,heap] = heap_pop(heap);
end
end

function heap = heap_initialize(capacity)
heap.key = zeros(capacity,1);
heap.task = zeros(capacity,1);
heap.size = 0;
end

function heap = heap_push(heap, key, task)
heap.size = heap.size + 1;
i = heap.size;
while i > 1
    parent = floor(i/2);
    if heap.key(parent) <= key
        break;
    end
    heap.key(i) = heap.key(parent);
    heap.task(i) = heap.task(parent);
    i = parent;
end
heap.key(i) = key;
heap.task(i) = task;
end

function [key, task, heap] = heap_pop(heap)
if heap.size == 0
    key = NaN;
    task = 0;
    return;
end
key = heap.key(1);
task = heap.task(1);
lastKey = heap.key(heap.size);
lastTask = heap.task(heap.size);
heap.size = heap.size - 1;
if heap.size == 0
    return;
end
i = 1;
while true
    left = 2*i;
    if left > heap.size
        break;
    end
    right = left+1;
    child = left;
    if right <= heap.size && heap.key(right) < heap.key(left)
        child = right;
    end
    if heap.key(child) >= lastKey
        break;
    end
    heap.key(i) = heap.key(child);
    heap.task(i) = heap.task(child);
    i = child;
end
heap.key(i) = lastKey;
heap.task(i) = lastTask;
end
