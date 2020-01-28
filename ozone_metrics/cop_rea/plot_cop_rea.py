# Clean up
plt.close('all')

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()


(data['O3'].sum(dim="time")/2).plot()

# Show it
plt.show(block=False)
