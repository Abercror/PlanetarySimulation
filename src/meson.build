src_files = [
    'main.cpp',
]

executable(
    'PlanetarySimulation',
    src_files,
    dependencies: [vulkan_dep, glfw_dep],
    include_directories: [glm_inc, include], 
    link_args : [
    '-framework', 'Cocoa',
    '-framework', 'IOKit',
    '-framework', 'CoreVideo'
    ],
    install : true
)