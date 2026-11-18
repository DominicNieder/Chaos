include("monodromy_method.jl")


# ==========================================================================================
#    perparing presentation: Understanding the study of phase space ->  surface of section
# ==========================================================================================
const SAVE_DATA_DIR_TRAJ   = joinpath(@__DIR__, "../../data/henon-heiles/presentation/poincare-map/traj/")
const SAVE_DATA_DIR_SEC   = joinpath(@__DIR__, "../../data/henon-heiles/presentation/poincare-map/map/")




# E= 0.1342
# trajectories chosen: chaotic, regular
E= 0.1342
p= (1.0,1.0,1.0)
sec_bound = boundary(E,p)
include("monodromy_method.jl")


const SAVE_DATA_DIR_TRAJ   = joinpath(@__DIR__, "../../data/henon-heiles/presentation/1trajectory/")
const SAVE_DATA_DIR_SEC   = joinpath(@__DIR__, "../../data/henon-heiles/presentation/poincare-map/map/")




# E= 0.1342
E= 0.1342
p= (1.0,1.0,1.0)
sec_bound = boundary(E,p)


# ===========================
#    first two trajectories chosen: chaotic, regular
# ===========================

# save data @
save_dir_2traj      = joinpath(@__DIR__, "../../data/henon-heiles/presentation/1trajectory/")
pathchaotic         = joinpath(save_dir_2traj, "chaotic_short")
pathreg             = joinpath(save_dir_2traj, "regular_short")

# chaotic
vchaotic = [-0.1, 0.0]
uchaotic = lift(vchaotic,E,p)
tchaotic = 500
# regular
vreg = [0.312, 0.0]
ureg = lift(vreg,E,p)
treg = 500

# trajectories
solchaotic = get_traj_and_map(uchaotic, tchaotic; dt=0.01)
solreg = get_traj_and_map(ureg, treg; dt=0.01)

#saving
save_trajectory_data(solchaotic, pathchaotic; format=:json)
save_trajectory_data(solreg, pathreg; format=:json)



# ======================================================
#    3D trajectories  -  complexity  -  to poincare map
# ======================================================

# Note: ti = 500 (if not chaotic && not poincaré section) 
# chaotic
v1 = [-0.1, 0.0]
u1 = lift(v1,E,p)
t1 = 100000
# regular
v2 = [0.312, 0.0]
u2 = lift(v2,E,p)
t2 = 500

# others for surface of section and showing complexity of trajectory space
v3 = [-0.025, 0.25]
u3 = lift(v3,E,p)
t3 = 500

v4 = [-0.13, 0.22]
u4 = lift(v4,E,p)
t4 = 500

v5 = [0.36, -0.16]
u5 = lift(v5,E,p)
t5 = 500

v6 = [0.1, 0.0]
u6 = lift(v6,E,p)
t6 = 500

v7 = [0.01, -0.25]
u7 = lift(v7,E,p)
t7 = 500

v8 = [-0.38, 0.1]
u8 = lift(v8,E,p)
t8 = 500


sol1 = get_traj_and_map(u1, t1; dt=0.01)
sol2 = get_traj_and_map(u2, t2; dt=0.01)
sol3 = get_traj_and_map(u3, t4; dt=0.01)
sol4 = get_traj_and_map(u4, t4; dt=0.01)
sol5 = get_traj_and_map(u5, t5; dt=0.01)
sol6 = get_traj_and_map(u6, t6; dt=0.01)
sol7 = get_traj_and_map(u7, t7; dt=0.01)
sol8 = get_traj_and_map(u8, t8; dt=0.01)


# chaecking section trajectory
py1  = sol1.sec_py
y1   = sol1.sec_y
py2  = sol2.sec_py
y2   = sol2.sec_y
py3  = sol3.sec_py
y3   = sol3.sec_y
py4  = sol4.sec_py
y4   = sol4.sec_y
py5  = sol5.sec_py
y5   = sol5.sec_y
py6  = sol6.sec_py
y6   = sol6.sec_y
py7  = sol7.sec_py
y7   = sol7.sec_y
py8  = sol8.sec_py
y8   = sol8.sec_y


set_style!(:dark)
fsec = Figure(size=(700,700))
axsec= Axis(fsec[1,1])
scatter!(axsec,sec_bound, color=INK)
scatter!(axsec, y1, py1, label="chaotic", color=pick_color(1))
scatter!(axsec, y2, py2, label="regular", color=pick_color(2))
scatter!(axsec, y3, py3, label="-", color=pick_color(3))
scatter!(axsec, y4, py4, label="-", color=pick_color(4))
scatter!(axsec, y5, py5, label="-", color=pick_color(5))
scatter!(axsec, y6, py6, label="-", color=pick_color(6))
scatter!(axsec, y7, py7, label="-", color=pick_color(7))
scatter!(axsec, y8, py8, label="-", color=pick_color(8))
display(fsec)


# trajectories
path1 = joinpath(SAVE_DATA_DIR_TRAJ, "1")
path2 = joinpath(SAVE_DATA_DIR_TRAJ, "2")
path3 = joinpath(SAVE_DATA_DIR_TRAJ, "3")
path4 = joinpath(SAVE_DATA_DIR_TRAJ, "4")
path5 = joinpath(SAVE_DATA_DIR_TRAJ, "5")
path6 = joinpath(SAVE_DATA_DIR_TRAJ, "6")
path7 = joinpath(SAVE_DATA_DIR_TRAJ, "7")
path8 = joinpath(SAVE_DATA_DIR_TRAJ, "8")

save_trajectory_data(sol1, path1; format=:json)
save_trajectory_data(sol2, path2; format=:json)
save_trajectory_data(sol3, path3; format=:json)
save_trajectory_data(sol4, path4; format=:json)
save_trajectory_data(sol5, path5; format=:json)
save_trajectory_data(sol6, path6; format=:json)
save_trajectory_data(sol7, path7; format=:json)
save_trajectory_data(sol8, path8; format=:json)

sec_path1 = joinpath(SAVE_DATA_DIR_SEC, "1")
sec_path2 = joinpath(SAVE_DATA_DIR_SEC, "2")
sec_path3 = joinpath(SAVE_DATA_DIR_SEC, "3")
sec_path4 = joinpath(SAVE_DATA_DIR_SEC, "4")
sec_path5 = joinpath(SAVE_DATA_DIR_SEC, "5")
sec_path6 = joinpath(SAVE_DATA_DIR_SEC, "6")
sec_path7 = joinpath(SAVE_DATA_DIR_SEC, "7")
sec_path8 = joinpath(SAVE_DATA_DIR_SEC, "8")

save_section_trajectory_data(sol1, sec_path1; format=:json)
save_section_trajectory_data(sol2, sec_path2; format=:json)
save_section_trajectory_data(sol3, sec_path3; format=:json)
save_section_trajectory_data(sol4, sec_path4; format=:json)
save_section_trajectory_data(sol5, sec_path5; format=:json)
save_section_trajectory_data(sol6, sec_path6; format=:json)
save_section_trajectory_data(sol7, sec_path7; format=:json)
save_section_trajectory_data(sol8, sec_path8; format=:json)



# one 3D trajectory, that x, y, px -> Pmap (show up to intersection 3) for all 

save_data_pmap_sec   = joinpath(@__DIR__, "../../data/henon-heiles/presentation/poincare-map/detail/sec/")
save_data_pmap_traj   = joinpath(@__DIR__, "../../data/henon-heiles/presentation/poincare-map/detail/map/")


tmax1 = sol1.sec_t[4]
tmax2 = sol2.sec_t[4]
tmax3 = sol3.sec_t[4]
tmax4 = sol4.sec_t[4]
tmax5 = sol5.sec_t[4]
tmax6 = sol6.sec_t[4]
tmax7 = sol7.sec_t[4]
tmax8 = sol8.sec_t[4]

get_int_time(t;fraction=2/7)  = t+ fraction*t

phi1 = get_traj_and_map(u1, get_int_time(tmax1); dt=0.01)
phi2 = get_traj_and_map(u2, get_int_time(tmax2); dt=0.01)
phi3 = get_traj_and_map(u3, get_int_time(tmax3); dt=0.01)
phi4 = get_traj_and_map(u4, get_int_time(tmax4); dt=0.01)
phi5 = get_traj_and_map(u5, get_int_time(tmax5); dt=0.01)
phi6 = get_traj_and_map(u6, get_int_time(tmax6); dt=0.01)
phi7 = get_traj_and_map(u7, get_int_time(tmax7); dt=0.01)
phi8 = get_traj_and_map(u8, get_int_time(tmax8); dt=0.01)



traj_path1 = joinpath(save_data_pmap_traj, "1")
traj_path2 = joinpath(save_data_pmap_traj, "2")
traj_path3 = joinpath(save_data_pmap_traj, "3")
traj_path4 = joinpath(save_data_pmap_traj, "4")
traj_path5 = joinpath(save_data_pmap_traj, "5")
traj_path6 = joinpath(save_data_pmap_traj, "6")
traj_path7 = joinpath(save_data_pmap_traj, "7")
traj_path8 = joinpath(save_data_pmap_traj, "8")

save_trajectory_data(phi1, traj_path1; format=:json)
save_trajectory_data(phi2, traj_path2; format=:json)
save_trajectory_data(phi3, traj_path3; format=:json)
save_trajectory_data(phi4, traj_path4; format=:json)
save_trajectory_data(phi5, traj_path5; format=:json)
save_trajectory_data(phi6, traj_path6; format=:json)
save_trajectory_data(phi7, traj_path7; format=:json)
save_trajectory_data(phi8, traj_path8; format=:json)

sec_path1 = joinpath(save_data_pmap_sec, "1")
sec_path2 = joinpath(save_data_pmap_sec, "2")
sec_path3 = joinpath(save_data_pmap_sec, "3")
sec_path4 = joinpath(save_data_pmap_sec, "4")
sec_path5 = joinpath(save_data_pmap_sec, "5")
sec_path6 = joinpath(save_data_pmap_sec, "6")
sec_path7 = joinpath(save_data_pmap_sec, "7")
sec_path8 = joinpath(save_data_pmap_sec, "8")

save_section_trajectory_data(sol1, sec_path1; format=:json)
save_section_trajectory_data(sol2, sec_path2; format=:json)
save_section_trajectory_data(sol3, sec_path3; format=:json)
save_section_trajectory_data(sol4, sec_path4; format=:json)
save_section_trajectory_data(sol5, sec_path5; format=:json)
save_section_trajectory_data(sol6, sec_path6; format=:json)
save_section_trajectory_data(sol7, sec_path7; format=:json)
save_section_trajectory_data(sol8, sec_path8; format=:json)




# ==================================================
#   surface of section: Scan through the energies 
# ==================================================


# energy range range(E0, 1/8, 1000)

# A, B, C trajectory at E0

# find bifurcating preriodic orbits for A, B, C theta=1/7, 1/6, 1/5, 1/4, 1/3,  

# find trajectories between A, B, C as funcitons of the poincare-map positions at E0
# i.e. range(abs(vA-vB),3), range(abs(vA-vC),3), range(abs(vB-vC),3)

# perform symmetry transformation to all periodic orbits




# ============================
#   finding periodic orbits
# ============================

# one 3D trajectory, that x, y, px -> Pmap