@testset "MAT output and input contracts" begin
    mktempdir() do dir
        mask=falses(3,2,2,4); mask[1,1,1,1]=mask[3,2,2,2]=true
        path=joinpath(dir,"regular.mat")
        tracks,info=track_mhw_4d_lonlat(mask;min_duration=1,save_file=path,add_idx5_alias=true)
        saved=matread(path)
        test_indexed_exact(indexed_from_mat(saved["tracks"],4),tracks)
        @test saved["info"]["n_tracks"]==info.n_tracks
        path=joinpath(dir,"empty.mat")
        track_mhw_4d_lonlat(falses(2,2,1,1);min_duration=1,save_file=path)
        @test haskey(matread(path),"tracks")
        matwrite(joinpath(dir,"mhw_ts.mat"),Dict("mhw_ts"=>reshape([NaN,0.,-1.,2.,Inf,3.,0.,1.],2,2,2)))
        matwrite(joinpath(dir,"lonlatoi.mat"),Dict("lon_used"=>[0.,180.],"lat_used"=>[-1.,1.]))
        d=load_test_data(dir;binary_rule=:positive)
        @test vec(d.mask)==[false,false,false,true,false,true,false,true]
        @test size(d.mask)==(2,2,2)
        @test count(load_test_data(dir;binary_rule=:finite).mask)==6
        @test isnan(load_test_data(dir).mask[1])
        # MATLAB scalar structs and cell arrays retain absolute days and order.
        original=Track3D([12329,12330],[[3.,1.],[2.,3.,1.]],[[1.,2.],[2.,1.,1.]],12329,3,[2],[12330])
        path=joinpath(dir,"overlap.mat")
        MHWTracking.save_tracks(path,[original],(n_tracks=1,))
        test_overlap_exact(load_matlab_tracks(path),[original])
        MHWTracking.save_tracks(path,[original,original],(n_tracks=2,))
        test_overlap_exact(load_matlab_tracks(path),[original,original])
        MHWTracking.save_tracks(path,Track3D[],(n_tracks=0,))
        @test isempty(load_matlab_tracks(path))
        matwrite(path,Dict("event"=>MHWTracking.matlab_struct(original)))
        test_overlap_exact(load_matlab_tracks(path;variable="event"),[original])
    end
end
