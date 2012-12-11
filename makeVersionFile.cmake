file(READ ${Kmernator_GIT_VERSION_FILE} Kmernator_VERSION_TAG_RAW LIMIT 1024)
string(STRIP "${Kmernator_VERSION_TAG_RAW}" Kmernator_VERSION_TAG)

set(Kmernator_VERSION "${Kmernator_VERSION_MAJOR}.${Kmernator_VERSION_MINOR}.${Kmernator_VERSION_TAG}")
configure_file(${Kmernator_VERSION_FILE_TEMPLATE} ${Kmernator_VERSION_FILE})
message("Building ${PROJECT_NAME} version ${Kmernator_VERSION}")
