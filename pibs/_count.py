import os


def countlines(start, lines=0, header=True, begin_start=None):
    if begin_start is None:
        begin_start = start

    if header:
        print("{:>10} |{:>10} | {:<20}".format("ADDED", "TOTAL", "FILE"))
        print("{:->11}|{:->11}|{:->20}".format("", "", ""))

    for thing in os.listdir(start):
        thing = os.path.join(start, thing)

        if os.path.isfile(thing):
            if thing.endswith(".py") and "_" not in thing:
                with open(thing, "r", encoding="utf-8") as f:
                    newlines = f.readlines()
                    newlines = len(newlines)
                    lines += newlines

                    reldir_of_thing = "." + thing.replace(begin_start, "")

                    print(
                        "{:>10} |{:>10} | {:<20}".format(
                            newlines, lines, reldir_of_thing
                        )
                    )

        elif os.path.isdir(thing):
            lines = countlines(thing, lines, header=False, begin_start=begin_start)

    # for thing in os.listdir(start):
    #     thing = os.path.join(start, thing)
    #     if os.path.isdir(thing):
    #         lines = countlines(thing, lines, header=False, begin_start=begin_start)

    return lines


if __name__ == "__main__":
    countlines(os.getcwd())
    input()
