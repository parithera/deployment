for db in "codeclarity" "knowledge" "config"; do
	docker compose -f docker-compose.yaml exec db sh -c "pg_dump -U postgres -d $db -Fc > ../../dump/$db.dump"
done