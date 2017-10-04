import numpy as np
from wordcloud import WordCloud

def plot_wordcloud(ax):

    baby = {"active": ["titanium dioxide", "zinc oxide"],
            "inactive": ["water", "capric triglyceride", "isodexadecane", "butyloctyl", "salicylate", "octyldodecyl citrate crosspolymer", "dimethicone", "sodium chloride", "ethylhexyl methoxycrylene", "caprylyl glycol", "PEG-8", "alumina", "glycerin", "sodium citrate", "tocopheryl acetate", "aloe barbadensis leaf extract", "phenoxyethanol"]}
    neutrogena = {"active": ["avobenzone", "homosalate", "octisalate", "octocrylene", "oxybenzone"],
                  "inactive": ["water", "butyloctyl salicylate", "styrene-acrylates copolymer", "silica", "diethylhexyl 2,6-napthalate", "vp-hexadecene copolymer", "dimethicone", "caprylyl methicone", "phenoxyethanol", "ethylhexyl glycerin", "glyceryl stearate", "PEG-100 stearate", "trimethylsiloxysilicate", "sodium polyacrylate", "methylparaden", "behenyl alcohol", "xanthan gum", "propylparaben", "ethylparaben", "fragrance", "disodium EDTA", "BHT", "iodopropynyl butylcarbamate"]}
    bullfrog = {"active": ["avobenzone", "homosalate", "octisalate", "octocrylene", "oxybenzone"],
                "inactive": ["acacia farnesiana extract", "acrylates-octylacrylamide copolymer", "aloe barbadensis leaf extract", "green tea", "chamomile flower extract", "cyclohexasiloxane", "cyclopentasiloxane", "ethylhexyl palmitate", "fragrance", "glycerin", "hydroxypropylcelluose", "isopropyl myristate", "lavender extract", "phenethyl benzoate", "PPG-12-SMDI coppolymer", "propylene glycol", "rosemary leaf extract", "SD alcohol 40", "tocopheryl acetate"]}
    badger = {"active": ["zinc oxide"],
              "inactive": ["sunflower seed oil", "beeswax", "jojoba seed oil", "tocopherol"]}
    coppertone = {"active": ["avobenzone", "homosalate", "octisalate", "octocrylene"],
                  "inactive": ["water", "aluminum starch octenylsuccinate", "styrene-acrylates copolymer", "glycerin", "polyester-27", "silica", "phenoxyethanol", "isododeance", "arachidyl alcohol", "beeswax", "ethylhexyl glycerin", "neopentyl", "glycol diheptanoate", "acrylates-C10-30 alkyl acrylate crosspolymer", "behenyl alcohol", "tocopherol", "arachidyl glucoside", "glyceryl stearate", "PEG-100 stearate", "potassium hydroxid", "disodium EDTA", "sdium ascorbyl phosphate", "fragrance"]}
    spray = {"active": ["avobenzone", "homosalate", "octisalate", "octocrylene"],
             "inactive": ["isobutane", "alcohol denat.", "isododecane", "diisopropyl adipate", "lauryl PEG-8 dimethicone", "phenylisopropyl dimethicone", "polyglycreyl-3 stearate", "isostearate-dimer dilinoleate copolymer", "capryl glycol", "methyl dihydroabietate", "fragrance", "ascorbyl palmitate", "tocopheryl acetate", "mineral oil", "panthenol", "water", "aloe barbadensis leaf extract"]}
    sunscreens = {'baby': baby,
                 'neutrogena': neutrogena,
                 'bullfrog': bullfrog,
                 'badger': badger,
                 'spray': spray,
                 'coppertone': coppertone}

    all_ingredients = []
    for sunscreen, ingredients in sunscreens.items():
        all_ingredients.extend(ingredients['active'])
        all_ingredients.extend(ingredients['inactive'])
    all_ingredients = np.unique(all_ingredients)

    fraction_stuck = {'baby': 1,
                      'neutrogena': 3,
                      'spray': 2,
                      'badger': -1,
                      'coppertone': -1,
                      'bullfrog': -1} # "danger"

    def get_ingredient_danger(ingredient):
        danger_points = []
        for sunscreen, ingredients in sunscreens.items():
            if ingredient in ingredients['active'] or ingredient in ingredients['inactive']:
                danger_points.append(fraction_stuck[sunscreen])
            else:
                for ix in np.hstack((ingredients['active'], ingredients['inactive'])):
                    if ingredient in ix.split(' '):
                        danger_points.append(fraction_stuck[sunscreen]*0.5)
        if len(danger_points) > 0:
            return np.sum(danger_points) #/ float(len(danger_points))
        else:
            return 0
            
    ingredient_dangers = {}
    for ingredient in all_ingredients:
        danger = get_ingredient_danger(ingredient)
        ingredient_dangers[ingredient] = danger

    ingr = ''
    for ingredient, danger in ingredient_dangers.items():
        if ingredient != 'water':
            danger += np.min(ingredient_dangers.values()) + 1
            i = ingredient
            i = i.replace(' ', '-')
            i = i.replace('-', '')
            i = i.replace(',', '')
            i = i.replace('.', '')
            ingr += (i + ' ')*int(danger)


    wordcloud = WordCloud(width=800*3, height=400*3, background_color="white").generate(ingr)

    ax.imshow(wordcloud)
    ax.axis("off")
