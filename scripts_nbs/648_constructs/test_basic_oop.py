import generate_648_cons as gen

PARTS_LINKERS_LIST = ["part1", "linker1", "part2", "linker2"]
basic_con = gen.BasicCon(*PARTS_LINKERS_LIST)
collection = gen.ConsCollection(basic_con)
