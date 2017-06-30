
from invite_filter import invite_filter

Genny = {
    "name":"Genny",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"no",
    "play tennis":"yes",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"yes",
    "exchange books":"super no",
    "play volleyball":"yes",
}

Christa = {
    "name":"Christa",
    "city":"Newport Beach",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"yes",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"no",
    "exchange books":"no",
    "play volleyball":"yes",
}

David = {
    "name":"David",
    "city":"Newport Beach",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"no",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"yes",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"yes",
    "exchange books":"yes",
    "play volleyball":"yes",
}

Sam_Delker = {
    "name":"Sam Delker",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"maybe",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"no",
    "exchange books":"no",
    "play volleyball":"maybe",
}

Sam_Little = {
    "name":"Sam Little",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"no",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"maybe",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"no",
    "exchange books":"maybe",
    "play volleyball":"maybe",
}

Miles = {
    "name":"Miles",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"maybe",
    "with someone who doesn't have an injured arm":"no",
    "is into engineering shit":"yes",
    "exchange books":"yes",
    "play volleyball":"yes",
}

Alex = {
    "name":"Alex",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"no",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"maybe",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"no",
    "exchange books":"yes",
    "play volleyball":"maybe",
}

Sam_Close = {
    "name":"Sam Close",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"yes",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"yes",
    "exchange books":"yes",
    "play volleyball":"yes",
}

Raphael = {
    "name":"Raphael",
    "city":"SB",
    "fighting":"no",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"no",
    "hang solo with someone":"no",
    "eat meat":"yes",
    "play tennis":"maybe",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"yes",
    "exchange books":"no",
    "play volleyball":"no",
}

Anil = {
    "name":"Anil",
    "city":"SB",
    "fighting":"yes",
    "bowl":"yes",
    "hike":"yes",
    "drink":"yes",
    "smoke weed":"yes",
    "hang solo with someone":"yes",
    "eat meat":"yes",
    "play tennis":"maybe",
    "with someone who doesn't have an injured arm":"yes",
    "is into engineering shit":"yes",
    "exchange books":"yes",
    "play volleyball":"yes",
}

#surf, beach day, hang out with a dog, book club

friends = [Genny, Christa, David, Sam_Delker, Sam_Little, Miles, Alex, Sam_Close, Raphael, Anil]



location = "SB"
list_of_what_you_want_to_do = []
list_of_what_you_want_to_do.append("smoke weed")
list_of_what_you_want_to_do.append("hike")


invite = invite_filter(friends,"city",location)
invite = invite_filter(invite,"fighting","no")

for what_you_want_to_do in list_of_what_you_want_to_do:
    invite = invite_filter(invite,what_you_want_to_do)
        
#print("If you want to ", list_of_what_you_want_to_do[0],end = "")
message = "If you want to "
#print("If you want to ")

count = 0
for what_you_want_to_do in list_of_what_you_want_to_do:
    #print(" and ",what_you_want_to_do,end="")
    if count == 0:
        message += what_you_want_to_do
    else:
        message += " and " + what_you_want_to_do
    count += 1

message += " in " + location + ", you can invite: \n"

print(message)
invite_names = []
for friend in invite:
    invite_names.append(friend['name'])

if invite_names:
    for name in invite_names:
        print(name)
else:
    print("Nobody. Sorry.")
