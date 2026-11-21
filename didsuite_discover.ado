*! didsuite_discover 0.1.0, 21nov2026, Matthew D. Webb; MIT license
program define didsuite_discover, rclass
    version 16
    local ids twfe did2s bjs flex csdid undid didint stacked lpdid
    local available
    foreach id in twfe did2s bjs flex csdid stacked lpdid {
        capture which fast`id'
        if !_rc {
            local available `available' `id'
            di as text "`id': fast`id' found; presence is not certification"
        }
        else di as text "`id': command missing"
    }
    foreach id in undid didint {
        capture which suite_`id'
        if !_rc local available `available' `id'
    }
    return local estimators "`ids'"
    return local available "`available'"
end

