project('PlantarySimulation', 'cpp', version: '0.1.0', default_options: ['cpp_std=c++17'])


vulkan_dep = dependency('vulkan', required : true)
glfw_dep = dependency('glfw3', required : true)

glm_inc = include_directories('/opt/homebrew/include')
include = include_directories('include')

subdir('src')
