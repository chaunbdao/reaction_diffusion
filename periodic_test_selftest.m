function periodic_test_selftest()

    assert(circular_height_for_test([1; 2; 9; 10], 10) == 4)
    assert(circular_height_for_test([3; 4; 5], 10) == 3)
    assert(circular_height_for_test(7, 10) == 1)

    [touches_initial,touches_final] = temporal_boundary_flags_for_test([2; 3; 4], 5);
    assert(~touches_initial)
    assert(~touches_final)

    [touches_initial,touches_final] = temporal_boundary_flags_for_test([1; 2; 3], 5);
    assert(touches_initial)
    assert(~touches_final)

    [touches_initial,touches_final] = temporal_boundary_flags_for_test([3; 4; 5], 5);
    assert(~touches_initial)
    assert(touches_final)

    space_time_area = [150; 150; 99; 150];
    touches_initial = [false; true; false; false];
    touches_final = [false; false; false; true];
    keep = keep_finished_patches_for_test(space_time_area, touches_initial, touches_final);
    assert(isequal(keep, [true; false; false; false]))

    disp('periodic_test_selftest passed')

end

function keep = keep_finished_patches_for_test(space_time_area, touches_initial, touches_final)
    keep = space_time_area >= 100 & ~touches_initial & ~touches_final;
end

function [touches_initial,touches_final] = temporal_boundary_flags_for_test(cols, runtime)
    touches_initial = any(cols == 1);
    touches_final = any(cols == runtime);
end

function height_pixels = circular_height_for_test(rows, points)
    rows = sort(rows(:));

    if isempty(rows)
        height_pixels = 0;
    elseif numel(rows) == 1
        height_pixels = 1;
    else
        gaps = diff(rows);
        wrap_gap = rows(1) + points - rows(end);
        height_pixels = points - max([gaps; wrap_gap]) + 1;
    end
end
