import generate_648_cons as gen
from dnabot import dnabot_app

PARTS_LINKERS_LIST = ["part1", "linker1", "part2", "linker2"]
basic_con = gen.BasicCon(*PARTS_LINKERS_LIST)
collection = gen.ConsCollection(basic_con)

constructs = gen.make_cons()
splitted_cons = constructs.split_cons(88)
print([len(build) for build in splitted_cons])

cons_lists = []
for build in splitted_cons:
    cons_lists.append([dnabot_app.process_construct(
        basic_con.parts_linkers) for basic_con in build])
print([dnabot_app.generate_clips_df(con_list) for con_list in cons_lists])