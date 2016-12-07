######################################################################
# Description
######################################################################

# 1. open Terminal
# 2. type R & hit enter
# 3. Paste and run (one line): 
  
nevek =c("Dóri", "Balázs", "Zsófi", "Marci", "Ebi", "Samu");
Addresses = c("Dóra Debreczeni <flydebrecenidora@gmail.com>", "Balázs Debreceni <dbalazs.gtd@gmail.com>", "Zsófi Debreceni <zsofiafilmes@gmail.com>", "Márton Debreceni <dmarci86@gmail.com>", "Ábel Vértesy <vertesy.abel@gmail.com>", "Sirajuddin (Sam) Cousins <soundsource@gmail.com>"); 
kap = sample(nevek); 

while (sum(kap == nevek)) {kap = sample(nevek)}; 
EMAILS = paste0("\n---------------\n\n",Addresses, "\n\nSzia ", nevek,"!\n\nAz idén te ", kap, " angyalkája leszel! \n\nJó karácsonyi készülődést: \nSecret Sánta" ); 
for(i in EMAILS) writeLines(i)


# nevek =c("Dóri", "Balázs", "Zsófi", "Marci", "Ebi", "Samu");Addresses = c("Dóra Debreczeni <flydebrecenidora@gmail.com>", "Balázs Debreceni <dbalazs.gtd@gmail.com>", "Zsófi Debreceni <zsofiafilmes@gmail.com>", "Márton Debreceni <dmarci86@gmail.com>", "Ábel Vértesy <vertesy.abel@gmail.com>", "Sirajuddin (Sam) Cousins <soundsource@gmail.com>"); kap = sample(nevek); while (sum(kap == nevek)) {kap = sample(nevek)}; EMAILS = paste0("\n---------------\n\n",Addresses, "\n\nSzia ", nevek,"!\n\nAz idén te ", kap, " angyalkája leszel! \n\nJó karácsonyi készülődést: \nSecret Sánta" ); for(i in EMAILS) writeLines(i)